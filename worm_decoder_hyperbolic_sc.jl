"""
worm_decoder_hyperbolic.jl
--------------------------
Implements the worm MCMC decoder for generic matchable qLDPC codes whose
decoding graph is specified directly as an adjacency list (rather than being
constructed from a parity check matrix).

For each syndrome realisation the decoder:
  1. Finds an initial matching using MWPM (via NetworkX).
  2. Runs worm MCMC sampling to estimate the probability of each logical coset.
  3. Returns the most probable logical correction.

Two decoders are benchmarked for each realisation:
  - MWPM                (initial matching, no MCMC)
  - Worm MLD            (MCMC majority-vote over sampled loop corrections)

Additionally, per-sample statistics are recorded:
  - Z_max    : fraction of samples in the most probable logical sector
  - Z_order  : fraction of samples in the correct logical sector

Simulation parameters (set in script):
  N_nodes    : number of decoding graph nodes
  adj        : adjacency list (1-indexed)
  logicals_list : list of logical operator edge sets (one per logical qubit)
  p_list     : list of bit-flip error rates
  N_real     : Monte Carlo realisations per p
  N_samp     : worm MCMC samples per syndrome
  t_therm    : thermalisation steps before sampling begins
  t_auto     : worm steps between samples (autocorrelation time)

Output:
  Results accumulated in p_MLD_list, p_MWPM_list, Z_max_list, Z_order_list
  per error rate p. Add CSV output as needed.

Dependencies:
  Random, LinearAlgebra, Graphs, PyCall (networkx)
"""

using Random
using LinearAlgebra
using Graphs
using PyCall

# =============================================================================
# Matching
# =============================================================================

# Find a minimum-weight perfect matching for the syndrome using MWPM (via
# NetworkX Blossom) and return the corresponding ±1 bond configuration.
#
# Note: any valid syndrome-consistent matching can be used as the worm starting
# point. A better initial match reduces thermalisation time but does not affect
# correctness.
function find_MWPM(synd, Adj_w, N_nodes, Adj)
    # Build weighted graph: edge weight = -log(odds ratio)
    Gr = nx.Graph()
    for (i, nbrs) in enumerate(Adj)
        for (k, j) in enumerate(nbrs)
            if j > i
                Gr.add_edge(i, j, weight=-log(Adj_w[i][k]))
            end
        end
    end

    # Build complete graph on syndrome nodes weighted by shortest-path distance
    CG = nx.Graph()
    for (idx, u) in enumerate(synd)
        for v in synd[idx+1:end]
            dist = nx.shortest_path_length(Gr, u, v, weight="weight")
            CG.add_edge(u, v, weight=-dist)
        end
    end

    # Run Blossom MWPM and recover paths in the original graph
    matching  = nx.max_weight_matching(CG, maxcardinality=true)
    path_list = []
    for (u, v) in matching
        path = nx.shortest_path(Gr, u, v, weight="weight")
        for i in 1:length(path)-1
            push!(path_list, [path[i], path[i+1]])
        end
    end

    # Convert traversed edges to a ±1 bond configuration
    matched_config = [fill(1, length(Adj[i])) for i in 1:N_nodes]
    for edge_ij in path_list
        node_i, node_j = edge_ij
        ind_ij = findfirst(x -> x == node_j, Adj[node_i])
        ind_ji = findfirst(x -> x == node_i, Adj[node_j])
        matched_config[node_i][ind_ij] = -1
        matched_config[node_j][ind_ji] = -1
    end
    return matched_config
end

# =============================================================================
# Graph Utilities (for initial matching fallback)
# =============================================================================

# Build a SimpleGraph from a dense adjacency matrix.
function build_graph_from_adjacency(A::Matrix{Int})
    g = SimpleGraph(size(A, 1))
    for i in 1:size(A, 1), j in i+1:size(A, 2)
        if A[i, j] != 0
            add_edge!(g, i, j)
        end
    end
    return g
end

# Find the shortest path between nodes a and b using Dijkstra's algorithm.
function find_path(A::Matrix{Int}, a::Int, b::Int)
    g  = build_graph_from_adjacency(A)
    if !has_path(g, a, b)
        return []
    end
    sp      = dijkstra_shortest_paths(g, a)
    path    = [b]
    current = b
    while current != a
        current = sp.parents[current]
        if current == 0
            return []
        end
        push!(path, current)
    end
    reverse!(path)
    return path
end

# Find a (random) perfect matching for syndrome_list by pairing nodes in
# shuffled order and routing each pair along a shortest path.
# Returns a ±1 bond configuration (mod-2 edge multiplicities).
function find_matching(N_nodes, A, syndrome_list)
    matching_matrix = zeros(Int, N_nodes, N_nodes)
    n = Int(length(syndrome_list) / 2)
    shuffle!(syndrome_list)
    for i in 1:n
        path_ab = find_path(A, syndrome_list[i], syndrome_list[n+i])
        for k in 1:(length(path_ab)-1)
            matching_matrix[path_ab[k], path_ab[k+1]] += 1
            matching_matrix[path_ab[k+1], path_ab[k]] += 1
        end
    end
    mm_mod_2      = matching_matrix .% 2
    matching_list = [mm_mod_2[i,:] for i in 1:N_nodes]
    return matching_list
end

# Convert a dense-matrix matching to a ±1 compact adjacency-list format
# compatible with the worm algorithm.
function generate_and_format_trial_error(N_nodes, A_matrix, syndrome)
    trial_matching = find_matching(N_nodes, A_matrix, syndrome)
    trial_matching_ = []
    for i in 1:N_nodes
        trial_i = [trial_matching[i][adj[i][j]] == 1 ? -1 : 1
                   for j in 1:length(adj[i])]
        push!(trial_matching_, trial_i)
    end
    return trial_matching_
end

# =============================================================================
# Error Generation and Syndrome Extraction
# =============================================================================

# Sample a bond disorder configuration at bit-flip rate p_FLIP.
# Returns A_quenched_dis where each entry is ±1, with -1 indicating a flipped
# bond. Each physical bond is sampled once and both directed edges are set
# consistently.
function generate_error(Adj_Mat, p_FLIP)
    A_quenched_dis = [fill(0, length(Adj_Mat[i])) for i in 1:length(Adj_Mat)]
    for i in 1:length(Adj_Mat)
        for j in 1:length(Adj_Mat[i])
            if A_quenched_dis[i][j] == 0
                val = rand() < p_FLIP ? -1 : 1
                A_quenched_dis[i][j] = val
                j_node    = Adj_Mat[i][j]
                i_index_j = findfirst(==(i), Adj_Mat[j_node])
                A_quenched_dis[j_node][i_index_j] = val
            end
        end
    end
    return A_quenched_dis
end

# Compute the error syndrome: returns the list of node indices where the
# product of incident bond variables is -1 (violated parity checks).
function find_syndrome(dis)
    synd_list = []
    for i in 1:length(dis)
        if prod(dis[i]) == -1
            push!(synd_list, i)
        end
    end
    return synd_list
end

# =============================================================================
# Worm Algorithm
# =============================================================================

# Advance the loop configuration by t_auto worm moves to produce a
# decorrelated sample.
function generate_independent_sample(original_loop_config, t_auto, A_weights,
                                      Adj_matrix, Quenched_Disorder)
    loop = copy(original_loop_config)
    for _ in 1:t_auto
        loop = generate_new_loop_config(loop, A_weights, Adj_matrix, Quenched_Disorder)
    end
    return loop
end

# Perform one worm move on the decoding graph.
#
# The worm starts at a random site, steps along bonds flipping each with
# Metropolis-Hastings probability, and closes when it returns to its start
# site — guaranteeing the updated config remains a valid loop satisfying all
# syndrome constraints.
#
# The degree-correction factor len(W[i])/len(W[j]) accounts for the
# non-uniform proposal distribution on graphs with varying node degrees.
function generate_new_loop_config(loop_config, Weights_Matrix, Adjacency_matrix,
                                   quenched_dis_matrix)
    N_sites      = length(Adjacency_matrix)
    start_site   = rand(1:N_sites)
    current_site = copy(start_site)

    while true
        direction = rand(1:length(Adjacency_matrix[current_site]))
        new_site  = Adjacency_matrix[current_site][direction]
        J_ij      = quenched_dis_matrix[current_site][direction]
        sigma_ij  = loop_config[current_site][direction]

        # Metropolis-Hastings acceptance weight with degree correction
        w_ij = if Int(0.5 * (1 + J_ij * sigma_ij)) == 0
            (1 / Weights_Matrix[current_site][direction]) *
            (length(Weights_Matrix[current_site]) / length(Weights_Matrix[new_site]))
        else
            Weights_Matrix[current_site][direction] *
            (length(Weights_Matrix[current_site]) / length(Weights_Matrix[new_site]))
        end
        w_ij = min(w_ij, 1.0)

        if rand() < w_ij
            loop_config[current_site][direction] = -sigma_ij
            old_index = findfirst(==(current_site), Adjacency_matrix[new_site])
            loop_config[new_site][old_index] = -sigma_ij
            current_site = new_site
        end

        current_site == start_site && break
    end

    return loop_config
end

# =============================================================================
# Topological Sector Measurement and MCMC Sampling
# =============================================================================

# Measure the logical sector of a loop configuration by evaluating each
# logical operator in logical_pairs_list.
# Returns a vector of ±1 values, one per logical qubit.
function topological_sector(loop_config, logical_pairs_list, adj)
    logical_state = []
    for logical_i in logical_pairs_list
        l_i_status = 1
        for l_i_edge in logical_i
            l_i_ind    = findfirst(==(l_i_edge[2]), adj[l_i_edge[1]])
            l_i_status *= loop_config[l_i_edge[1]][l_i_ind]
        end
        push!(logical_state, l_i_status)
    end
    return logical_state
end

# Run the worm MCMC sampler. Thermalises for t_THERMALIZE steps, then draws
# N_samples decorrelated loop configurations separated by t_AUTO steps.
# Tallies the number of times each distinct logical sector is visited.
#
# Returns (logical_states_list, logical_states_tally), where each entry of
# logical_states_list is a distinct sector vector and the corresponding tally
# is the number of samples in that sector. The most probable sector is
# argmax(tally).
function run_sampling(initial_loop, N_samples, t_THERMALIZE, t_AUTO, A_w, A_,
                      Quenched_Disorder, logical_pairs_list)
    loop = generate_independent_sample(initial_loop, t_THERMALIZE, A_w, A_, Quenched_Disorder)
    logical_states_list  = []
    logical_states_tally = []

    for _ in 1:N_samples
        loop           = generate_independent_sample(loop, t_AUTO, A_w, A_, Quenched_Disorder)
        logical_states = topological_sector(loop, logical_pairs_list, A_)

        if logical_states in logical_states_list
            logical_states_tally[findfirst(==(logical_states), logical_states_list)] += 1
        else
            push!(logical_states_list,  logical_states)
            push!(logical_states_tally, 1)
        end
    end

    return logical_states_list, logical_states_tally
end

# =============================================================================
# Graph Utilities
# =============================================================================

# Build a compact adjacency list from a list of directed edge pairs.
# Edge indices are shifted from 0-based to 1-based.
function gen_original_lattice_adj(edges, num_v)
    adj = [Int64[] for _ in 1:num_v]
    for e_pair in edges
        push!(adj[e_pair[1]+1], e_pair[2]+1)
        push!(adj[e_pair[2]+1], e_pair[1]+1)
    end
    return adj
end

# Compute the symmetric difference (element-wise product) of two ±1 bond
# configurations. A bond is -1 iff exactly one of the two configs has it as -1,
# corresponding to composition of two Pauli error chains.
function symmetric_dif(a1, a2)
    return [[a1[j][i] * a2[j][i] for i in 1:length(a1[j])] for j in 1:length(a1)]
end

# =============================================================================
# Simulation Parameters
# =============================================================================

t_therm = 10000     # thermalisation steps before sampling
t_auto  = 1000      # worm steps between samples (autocorrelation time)
N_real  = 50        # Monte Carlo realisations per p
N_samp  = 10000     # worm MCMC samples per syndrome

p_list = [0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025,
          0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045]

# =============================================================================
# Decoding Graph and Logical Operators
# Full data files are in data/hyperbolic_decoding_graphs/
# =============================================================================

# adj = [...]           # adjacency list
# loops_list = [...]    # logical operator edge sets

# =============================================================================
# Main Simulation Loop
# =============================================================================
N_nodes = length(adj)
nx           = pyimport("networkx")
k_logicals   = length(logicals_list)
trivial_state = fill(1, k_logicals)

original_loop_config = [fill(1, length(adj[i])) for i in 1:N_nodes]

for p in p_list
    a_weights = [fill(p / (1 - p), length(adj[i])) for i in 1:N_nodes]

    p_MLD_list   = []
    p_MWPM_list  = []
    Z_max_list   = []
    Z_order_list = []

    for _ in 1:N_real
        error_hyp = generate_error(adj, p)
        synd      = find_syndrome(error_hyp)

        # --- Step 1: MWPM initial matching ---
        MWPM_sol = find_MWPM(synd, a_weights, N_nodes, adj)

        # --- Step 2: Worm MCMC sampling ---
        results = run_sampling(original_loop_config, N_samp, t_therm, t_auto,
                               a_weights, adj, MWPM_sol, logicals_list)
        j_max      = argmax(results[2])
        correction = results[1][j_max]

        push!(Z_max_list, results[2][j_max] / N_samp)

        # --- Step 3: Evaluate correctness ---
        loop_c         = symmetric_dif(error_hyp, MWPM_sol)
        correct_logical = topological_sector(loop_c, logicals_list, adj)

        push!(p_MLD_list,  correction == correct_logical ? 1 : 0)
        push!(p_MWPM_list, correct_logical == trivial_state ? 1 : 0)

        if correct_logical in results[1]
            correct_log_index = findfirst(==(correct_logical), results[1])
            push!(Z_order_list, results[2][correct_log_index] / N_samp)
        else
            push!(Z_order_list, 0)
        end
    end
end
