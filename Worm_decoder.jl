"""
worm_decoder.jl
---------------
Implements the Worm MCMC decoder for matchable qLDPC codes.

Given a syndrome (set of violated detectors) and a decoding graph, the decoder:
  1. Finds an initial configuration using MWPM (via NetworkX).
  2. Estimates the probability of each logical coset by sampling.
  3. Returns the most likely logical correction and its estimated probability.

Inputs (command-line arguments):
    ARGS[1] : path to DEM JSON file   (adjacency list, logical edges, MCMC parameters)
    ARGS[2] : path to syndrome JSON file (list of syndromes to decode)

Output (stdout):
    Two Julia Any[] arrays on separate lines:
        Line 1 : per-syndrome correction decisions (±1)
        Line 2 : per-syndrome MLD probability estimates

Dependencies:
    JSON, PyCall (networkx), Random, LinearAlgebra, DataFrames, Graphs
"""

using JSON
using Random
using LinearAlgebra
using DataFrames
using Graphs
using PyCall


# =============================================================================
# Worm Algorithm
# =============================================================================

"""
    generate_new_loop_config(loop_config, Weights_Matrix, Adjacency_matrix,
                             quenched_dis_matrix)

Perform one worm move: propose a bond flip along a stochastic walk starting from a
random site, accepting each step with a Metropolis-Hastings weight.

The walk terminates and the loop closes when the worm returns to its starting
site.

# Arguments
- `loop_config`         : current bond configuration (list of lists of ±1)
- `Weights_Matrix`      : bond weights (log-odds of error on each bond)
- `Adjacency_matrix`    : adjacency list (1-indexed)
- `quenched_dis_matrix` : quenched disorder (reference matching)

# Returns
- Updated `loop_config`
"""
function generate_new_loop_config(loop_config, Weights_Matrix, Adjacency_matrix, quenched_dis_matrix)
    new_loop_config = false
    N_sites = length(Adjacency_matrix)
    start_site = rand(1:N_sites)
    current_site = copy(start_site)

    while !new_loop_config
        direction = rand(1:length(Adjacency_matrix[current_site]))
        new_site  = Adjacency_matrix[current_site][direction]
        J_ij      = quenched_dis_matrix[current_site][direction]
        sigma_ij  = loop_config[current_site][direction]

        # Metropolis-Hastings acceptance weight
        if Int(0.5 * (1 + J_ij * sigma_ij)) == 0
            w_ij = (1 / Weights_Matrix[current_site][direction]) *
                   (length(Weights_Matrix[current_site]) / length(Weights_Matrix[new_site]))
        else
            w_ij = Weights_Matrix[current_site][direction] *
                   (length(Weights_Matrix[current_site]) / length(Weights_Matrix[new_site]))
        end
        w_ij = min(w_ij, 1.0)

        # Accept or reject the proposed bond flip
        if rand() < w_ij
            loop_config[current_site][direction] = -sigma_ij
            old_index = findfirst(==(current_site), Adjacency_matrix[new_site])
            loop_config[new_site][old_index] = -sigma_ij
            current_site = new_site
        end

        # Loop closes when worm returns to start
        if current_site == start_site
            new_loop_config = true
        end
    end

    return loop_config
end


"""
    generate_independent_sample(original_loop_config, t_auto, A_weights,
                                Adj_matrix, Quenched_Disorder)

Advance the loop configuration by `t_auto` worm moves to produce an
approximately independent sample.

# Returns
- Updated loop configuration after `t_auto` steps.
"""
function generate_independent_sample(original_loop_config, t_auto, A_weights, Adj_matrix, Quenched_Disorder)
    loop = copy(original_loop_config)
    for _ in 1:t_auto
        loop = generate_new_loop_config(loop, A_weights, Adj_matrix, Quenched_Disorder)
    end
    return loop
end


# =============================================================================
# Topological Sector
# =============================================================================

"""
    topological_sector(loop_config, logical_pairs_list, adj)

Compute the topological (logical) sector of a loop configuration by evaluating
the product of bond variables along each logical edge.

Currently supports a single logical observable (k=1).

# Returns
- `logical_state` : list of ±1, one entry per logical observable.
"""
function topological_sector(loop_config, logical_pairs_list, adj)
    k = 1   # number of logical observables (generalise for k > 1 if needed)
    logical_state = []

    for i in 1:k
        l_i_status = 1
        for l_i_edge in logical_pairs_list
            l_i_ind    = findfirst(==(l_i_edge[2]), adj[l_i_edge[1]])
            l_i_status = l_i_status * loop_config[l_i_edge[1]][l_i_ind]
        end
        push!(logical_state, l_i_status)
    end

    return logical_state
end


# =============================================================================
# MCMC Sampling
# =============================================================================

"""
    run_sampling(initial_loop, N_samples, t_THERMALIZE, t_AUTO,
                 A_w, A_, Quenched_Disorder, logical_pairs_list)

Thermalise the loop configuration, then draw `N_samples` independent samples
and count the number of times the loop is in each logical sector.

# Returns
- `logical_state_counter` : [count_sector_+1, count_sector_-1]
"""
function run_sampling(initial_loop, N_samples, t_THERMALIZE, t_AUTO, A_w, A_, Quenched_Disorder, logical_pairs_list)
    # Thermalise before collecting samples
    loop = generate_independent_sample(initial_loop, t_THERMALIZE, A_w, A_, Quenched_Disorder)

    logical_state_counter = [0, 0]
    for _ in 1:N_samples
        loop          = generate_independent_sample(loop, t_AUTO, A_w, A_, Quenched_Disorder)
        logical_states = topological_sector(loop, logical_pairs_list, A_)

        # Tally logical sector (generalise for multiple logicals if needed)
        if logical_states[1] == 1
            logical_state_counter[1] += 1
        else
            logical_state_counter[2] += 1
        end
    end

    return logical_state_counter
end


# =============================================================================
# Graph Utilities
# =============================================================================

"""
    find_MWPM(synd, Adj_w, N_nodes, Adj)

Find a perfect matching for a given syndrome using MWPM (via NetworkX) and
return the corresponding bond configuration.

Note: any valid matching can be used here — the worm algorithm only requires a
syndrome-consistent starting configuration. A better initial match may reduce
thermalisation time, but does not affect correctness.

# Arguments
- `synd`    : list of syndrome node indices (1-indexed)
- `Adj_w`   : bond weight adjacency list (log-odds)
- `N_nodes` : total number of nodes
- `Adj`     : adjacency list (1-indexed)

# Returns
- `matched_config` : bond configuration (list of lists of ±1)
"""
function find_MWPM(synd, Adj_w, N_nodes, Adj)
    # Build weighted graph in NetworkX
    Gr = nx.Graph()
    for (i, nbrs) in enumerate(Adj)
        for (k, j) in enumerate(nbrs)
            if j > i    # avoid adding each edge twice
                w = -log(Adj_w[i][k])
                Gr.add_edge(i, j, weight=w)
            end
        end
    end

    # Build complete graph on syndrome nodes with shortest-path distances
    CG = nx.Graph()
    for (idx, u) in enumerate(synd)
        for v in synd[idx+1:end]
            dist = nx.shortest_path_length(Gr, u, v, weight="weight")
            CG.add_edge(u, v, weight=-dist)
        end
    end

    # Find minimum-weight perfect matching
    matching = nx.max_weight_matching(CG, maxcardinality=true)

    # Recover the matched bond configuration from the matching paths
    path_list = []
    for (u, v) in matching
        path = nx.shortest_path(Gr, u, v, weight="weight")
        for i in 1:length(path)-1
            push!(path_list, [path[i], path[i+1]])
        end
    end

    matched_config = [fill(1, length(Adj[i])) for i in 1:N_nodes]
    for edge_ij in path_list
        node_i = edge_ij[1]
        node_j = edge_ij[2]
        ind_ij = findfirst(x -> x == node_j, Adj[node_i])
        ind_ji = findfirst(x -> x == node_i, Adj[node_j])
        matched_config[node_i][ind_ij] = -1
        matched_config[node_j][ind_ji] = -1
    end

    return matched_config
end


# =============================================================================
# Main: Load Inputs and Run Decoder
# =============================================================================

# --- Parse command-line arguments ---
filename_DEM  = ARGS[1]
filename_synd = ARGS[2]

# --- Load DEM and syndrome data from JSON ---
data_DEM  = JSON.parsefile(filename_DEM)
data_synd = JSON.parsefile(filename_synd)

# --- Extract DEM fields ---
logicals     = data_DEM["logicals"]
adj          = data_DEM["adj"]
adj_w        = data_DEM["adj_w"]
t_AUTO       = data_DEM["t_autocorrelation"]
t_THERMALIZE = data_DEM["t_thermalize"]
N_samples    = data_DEM["N_samples"]

# --- Extract syndromes and delete the temporary syndrome file ---
syndrome_list = data_synd["syndromes"]
rm(filename_synd; force=true)

# --- Initialise NetworkX (called via PyCall for MWPM) ---
nx = pyimport("networkx")

N_nodes = length(adj)

correction_output_list = []
correction_prob_list   = []

# =============================================================================
# Decode each syndrome
# =============================================================================

for syndrome in syndrome_list

    # 1. Find an initial loop configuration satisfying the syndrome constraints.
    #    MWPM is used here, but any valid matching can be substituted —
    #    the worm algorithm does not depend on the quality of the initial match.
    reference_matching_ = find_MWPM(syndrome, adj_w, N_nodes, adj)

    # 2. Determine the logical sector of the initial MWPM solution
    initial_log = topological_sector(reference_matching_, logicals, adj)

    # 3. Start the worm MCMC from the trivial (all +1) loop configuration
    initial_loop = [fill(1, length(adj[i])) for i in 1:N_nodes]

    # 4. Run MCMC sampling
    logical_sampling_results = run_sampling(
        initial_loop, N_samples, t_THERMALIZE, t_AUTO,
        adj_w, adj, reference_matching_, logicals
    )

    # 5. Pick the most probable logical sector and record the correction
    if logical_sampling_results[1] > N_samples / 2
        l = 1
        push!(correction_prob_list, logical_sampling_results[1] / N_samples)
    else
        l = -1
        push!(correction_prob_list, logical_sampling_results[2] / N_samples)
    end

    push!(correction_output_list, l * initial_log[1])
end

# --- Output results to stdout (read by the Python driver script) ---
println(correction_output_list)
println(correction_prob_list)
