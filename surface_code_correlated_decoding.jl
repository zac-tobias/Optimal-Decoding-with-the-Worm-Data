using Random
using LinearAlgebra
using CSV
using DataFrames
using Graphs
using PyCall

# =============================================================================
# RSC Code Construction
# =============================================================================

# Return the X stabiliser check supports and the X logical operator qubit list
# for an L×L rotated surface code (RSC). L must be odd.
#
# The final entry of x_checks_list is a virtual boundary node whose support
# spans all left- and right-boundary qubits, allowing the matching to pair
# syndrome nodes near opposite code boundaries without explicit correction.
function find_X_checks(L)
    x_checks_list = []
    x_logical = []

    # Bulk X checks (weight-4 plaquettes)
    for i in 1:(L-1)/2
        x = Int(2*i - 1)
        for j in 1:(L-1)/2
            y = Int(2*j - 1)
            push!(x_checks_list, [coords_to_qubit(x+1,y,L), coords_to_qubit(x+2,y,L),
                                   coords_to_qubit(x+1,y+1,L), coords_to_qubit(x+2,y+1,L)])
            push!(x_checks_list, [coords_to_qubit(x,y+1,L), coords_to_qubit(x,y+2,L),
                                   coords_to_qubit(x+1,y+1,L), coords_to_qubit(x+1,y+2,L)])
        end
    end

    # Boundary X checks (weight-2)
    for k in 1:(L-1)/2
        x = Int(2*k - 1)
        push!(x_checks_list, [coords_to_qubit(x,1,L),   coords_to_qubit(x+1,1,L)])
        push!(x_checks_list, [coords_to_qubit(x+1,L,L), coords_to_qubit(x+2,L,L)])
    end

    # Virtual boundary node: all left- and right-boundary qubits
    boundary_node_x = []
    for y in 1:L
        push!(boundary_node_x, coords_to_qubit(1,y,L))
        push!(x_logical,       coords_to_qubit(1,y,L))
        push!(boundary_node_x, coords_to_qubit(L,y,L))
    end
    push!(x_checks_list, boundary_node_x)

    return x_checks_list, x_logical
end

# Return the Z stabiliser check supports for an L×L RSC. L must be odd.
# The final entry is a virtual boundary node spanning top and bottom boundaries.
function find_Z_checks(L)
    z_checks_list = []

    # Bulk Z checks (weight-4 plaquettes)
    for i in 1:(L-1)/2
        x = Int(2*i - 1)
        for j in 1:(L-1)/2
            y = Int(2*j - 1)
            push!(z_checks_list, [coords_to_qubit(x,y,L),   coords_to_qubit(x+1,y,L),
                                   coords_to_qubit(x,y+1,L), coords_to_qubit(x+1,y+1,L)])
            push!(z_checks_list, [coords_to_qubit(x+1,y+1,L), coords_to_qubit(x+1,y+2,L),
                                   coords_to_qubit(x+2,y+1,L), coords_to_qubit(x+2,y+2,L)])
        end
    end

    # Boundary Z checks (weight-2)
    for k in 1:(L-1)/2
        y = Int(2*k - 1)
        push!(z_checks_list, [coords_to_qubit(L,y,L),   coords_to_qubit(L,y+1,L)])
        push!(z_checks_list, [coords_to_qubit(1,y+1,L), coords_to_qubit(1,y+2,L)])
    end

    # Virtual boundary node: all top- and bottom-boundary qubits
    boundary_node_z = []
    for x in 1:L
        push!(boundary_node_z, coords_to_qubit(x,1,L))
        push!(boundary_node_z, coords_to_qubit(x,L,L))
    end
    push!(z_checks_list, boundary_node_z)

    return z_checks_list
end

# Map 2D lattice coordinates (x, y) to a qubit index (row-major, 1-indexed).
function coords_to_qubit(x, y, L)
    return x + (y-1)*L
end

# Build the binary parity check matrix H ∈ {0,1}^{m×n} from check supports.
# H[i,q] = 1 iff qubit q is in the support of check i.
function find_check_mat(checks_list, n)
    m = length(checks_list)
    H = zeros(Int, m, n)
    for i in 1:m
        for qu in checks_list[i]
            H[i, qu] = 1
        end
    end
    return H
end

# =============================================================================
# Decoding Graph Construction
# =============================================================================

# Convert a dense adjacency matrix to a compact adjacency list.
# Multi-edges (adj[i,k] > 1) are preserved by repeating the neighbour index,
# arising when two stabilisers share more than one qubit in their supports.
function adj_to_compact_adj(adj)
    N = size(adj, 1)
    adj_compact = []
    for i in 1:N
        adj_i_compact = []
        for k in 1:N
            if adj[i,k] > 0 && k != i
                for _ in 1:adj[i,k]
                    push!(adj_i_compact, k)
                end
            end
        end
        push!(adj_compact, adj_i_compact)
    end
    return adj_compact
end

# Build the edge→qubit map for the decoding graph.
# Each edge (m1, m2) corresponds to the qubit in the intersection of the
# supports of checks m1 and m2 (i.e. H[m1,q] = H[m2,q] = 1).
# Returns qubit_list_mapping[m1][k] = qubit for the k-th neighbour of m1.
function edge_to_qubit(adj_comp, H, n)
    qubit_list_mapping = []
    N   = length(adj_comp)
    H_q = transpose(H)
    for m_1 in 1:N
        qubit_list_i = []
        for m_2 in adj_comp[m_1]
            for i in 1:n
                h_q_i_ind = findall(x -> x == 1, H_q[i,:])
                if m_1 in h_q_i_ind && m_2 in h_q_i_ind && i ∉ qubit_list_i
                    push!(qubit_list_i, i)
                end
            end
        end
        push!(qubit_list_mapping, qubit_list_i)
    end
    return qubit_list_mapping
end

# Build the inverse qubit→edge map.
# qubit_edges[q] is a list of [node_index, edge_index] pairs identifying every
# decoding graph edge that corresponds to qubit q. Used during worm updates to
# locate and flip the reverse edge when the worm crosses a bond.
function qubits_to_edge(qubit_mapping, n_qubits)
    qubit_edges = [[] for _ in 1:n_qubits]
    for i in 1:length(qubit_mapping)
        for j in 1:length(qubit_mapping[i])
            push!(qubit_edges[qubit_mapping[i][j]], [i, j])
        end
    end
    return qubit_edges
end

# =============================================================================
# Error Generation and Syndrome Extraction
# =============================================================================

# Sample an i.i.d. depolarizing error on n qubits at physical rate p.
# Returns binary qubit_list_x and qubit_list_z (1 = error present).
function generate_error(p, n)
    qubit_list_x = Int[]
    qubit_list_z = Int[]
    for _ in 1:n
        r = rand()
        if r < p/3           # Y error: both X and Z components
            push!(qubit_list_x, 1); push!(qubit_list_z, 1)
        elseif r < 2p/3      # X error
            push!(qubit_list_x, 1); push!(qubit_list_z, 0)
        elseif r < p         # Z error
            push!(qubit_list_x, 0); push!(qubit_list_z, 1)
        else                  # No error
            push!(qubit_list_x, 0); push!(qubit_list_z, 0)
        end
    end
    return qubit_list_x, qubit_list_z
end

# Compute the error syndrome from a qubit error pattern.
# Returns the list of check indices where the parity of errors in the support
# is odd (i.e. violated stabilisers).
function find_syndrome(qubit_mapping, qubits)
    syndrome = []
    for i in 1:length(qubit_mapping)
        parity = 1
        for q in qubit_mapping[i]
            if qubits[q] == 1
                parity = -parity
            end
        end
        if parity == -1
            push!(syndrome, i)
        end
    end
    return syndrome
end

# Convert a qubit-level X error pattern to a ±1 bond configuration on the
# X decoding graph. Each errored qubit flips the two graph edges it connects.
function x_error_edges(q_x_errors, adj_x, q_edge_map)
    N_x = length(adj_x)
    x_error_edge_config = [fill(1, length(adj_x[i])) for i in 1:N_x]
    for i in 1:length(q_x_errors)
        if q_x_errors[i] == 1
            q_edges = q_edge_map[i]
            x_error_edge_config[q_edges[1][1]][q_edges[1][2]] = -1
            x_error_edge_config[q_edges[2][1]][q_edges[2][2]] = -1
        end
    end
    return x_error_edge_config
end

# =============================================================================
# Logical Operator Measurement
# =============================================================================

# Find the edge indices within the virtual boundary node that correspond to
# X logical operator qubits. Used by logical_measure to evaluate the logical
# outcome of a correction.
function find_X_logicals(x_log_qubit, edge_to_qubits)
    x_log_ind = []
    boundary_qubits = edge_to_qubits[end]
    for q in x_log_qubit
        if q in boundary_qubits
            push!(x_log_ind, findfirst(==(q), boundary_qubits))
        end
    end
    return x_log_ind
end

# Evaluate the X logical operator by multiplying bond variables along the
# logical cycle within the virtual boundary node.
# Returns +1 (no logical error) or -1 (logical error).
function logical_measure(X_config, log_indices)
    l = 1
    for i in log_indices
        l *= X_config[end][i]
    end
    return l
end

# =============================================================================
# Matching
# =============================================================================

# Find a minimum-weight perfect matching for the syndrome using MWPM (via
# NetworkX Blossom) and return the corresponding ±1 bond configuration.
#
# Note: any valid syndrome-consistent matching can be used as the worm starting
# point. A better initial match reduces thermalisation time but does not affect
# correctness.
function find_MWPM(synd, Adj_w, N_nodes, adj)
    # Build weighted graph: edge weight = -log(odds ratio)
    Gr = nx.Graph()
    for (i, nbrs) in enumerate(adj)
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
    matched_config = [fill(1, length(adj[i])) for i in 1:N_nodes]
    for edge_ij in path_list
        node_i, node_j = edge_ij
        ind_ij = findfirst(x -> x == node_j, adj[node_i])
        ind_ji = findfirst(x -> x == node_i, adj[node_j])
        matched_config[node_i][ind_ij] = -1
        matched_config[node_j][ind_ji] = -1
    end
    return matched_config
end

# =============================================================================
# Worm Algorithm
# =============================================================================

# Advance the loop configuration by t_auto worm moves to produce a
# decorrelated sample.
function generate_independent_sample(original_loop_config, t_auto, A_weights,
                                      Adj_matrix, Quenched_Disorder,
                                      qubit_mapping, qubit_edges)
    loop = deepcopy(original_loop_config)
    for _ in 1:t_auto
        loop = generate_new_loop_config(loop, A_weights, Adj_matrix,
                                        Quenched_Disorder, qubit_mapping, qubit_edges)
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
                                   quenched_dis_matrix, edge_to_qubit, qubit_to_edge)
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

            # Locate the reverse edge at new_site via qubit_to_edge and flip it
            qubit_edges_ij = qubit_to_edge[edge_to_qubit[current_site][direction]]
            old_index = qubit_edges_ij[1][1] == new_site ? qubit_edges_ij[1][2] : qubit_edges_ij[2][2]
            loop_config[new_site][old_index] = -sigma_ij

            current_site = new_site
        end

        current_site == start_site && break
    end

    return loop_config
end

# Run the worm MCMC decoder. Draws n_samples loop configurations, tallies the
# logical sector of (error ∘ matching ∘ worm_loop) for each, and accumulates
# per-edge marginal error probabilities for use in the next parsing iteration.
#
# Returns (success::Int, error_dist) where success = 1 if the trivial sector
# is most probable (ties broken uniformly at random), and error_dist[i][j] is
# the estimated probability that edge (i,j) carries an error.
function worm_decoder(adj_w, adj, matching, n_samples, t_auto,
                      edges_to_qubits, qubit_edges, error_chain, log_ind)
    sector_tally  = [0, 0]
    matching_loop = [fill(1, length(adj[i])) for i in 1:length(adj)]
    error_dist    = [fill(0, length(adj[i])) for i in 1:length(adj)]
    N_nodes       = length(adj)

    for _ in 1:n_samples
        matching_loop = generate_independent_sample(
            matching_loop, t_auto, adj_w, adj, deepcopy(matching),
            edges_to_qubits, qubit_edges
        )
        matching_MLD = symmetric_dif(matching, matching_loop)

        # Accumulate marginal error probability: p(edge error) = (1 - σ) / 2
        error_dist = [(error_dist[i] .+ (-matching_MLD[i] .+ 1) / (2*n_samples))
                      for i in 1:N_nodes]

        total_op    = symmetric_dif(error_chain, matching_MLD)
        logical_val = logical_measure(total_op, log_ind)
        logical_val == 1 ? sector_tally[1] += 1 : sector_tally[2] += 1
    end

    result = if sector_tally[1] == sector_tally[2]
        rand([0, 1])
    elseif sector_tally[1] > sector_tally[2]
        1
    else
        0
    end

    return result, error_dist
end

# =============================================================================
# Graph Utilities
# =============================================================================

# Compute the symmetric difference (element-wise product) of two ±1 bond
# configurations. A bond is -1 iff exactly one of the two configs has it as -1,
# corresponding to composition of two Pauli error chains.
function symmetric_dif(a1, a2)
    return [[a1[j][i] * a2[j][i] for i in 1:length(a1[j])] for j in 1:length(a1)]
end

# Find all nodes where the product of incident bond variables is -1.
# Used to verify that a matching correctly resolves a syndrome.
function find_syndrome_edge_conf(dis)
    synd_list = []
    for i in 1:length(dis)
        if prod(dis[i]) == -1
            push!(synd_list, i)
        end
    end
    return synd_list
end

# =============================================================================
# Correlated Reweighting
# =============================================================================

# Reweight the decoding graph for one Pauli type using soft marginal error
# probabilities from the worm decoder for the other Pauli type.
#
# This accounts for Y-error correlations between X and Z syndromes under
# depolarizing noise. For an edge crossing qubit q with estimated error
# probability p_q (from the opposite decoder):
#   p_ij = 0.5·p_q + (1 - p_q)·(p/3)/(1 - 2p/3)
# The first term is the Y-error contribution (correlated X and Z);
# the second is the independent contribution conditional on no Y error.
#
# Arguments:
#   edge_prob_list  — per-edge marginal error probs from the source graph
#   qubit_mapping_1 — edge→qubit map for the source graph
#   qubit_mapping_2 — edge→qubit map for the target graph to be reweighted
#   p               — marginal single-qubit X (or Z) rate = 2·p_depol/3
function reweight_XZ(edge_prob_list, qubit_mapping_1, qubit_mapping_2, L, p, N_1, N_2)
    # Accumulate per-qubit error probabilities from the source graph
    qubit_probs = zeros(L^2)
    for i in 1:N_1
        for j in 1:length(edge_prob_list[i])
            qubit_probs[qubit_mapping_1[i][j]] = edge_prob_list[i][j]
        end
    end

    # Reweight each edge of the target graph using its qubit's marginal prob
    adj_w_2 = []
    for i in 1:N_2
        adj_w_i = []
        for j in 1:length(qubit_mapping_2[i])
            p_q  = qubit_probs[qubit_mapping_2[i][j]]
            p_ij = 0.5*p_q + (1 - p_q)*(p/3)/(1 - 2p/3)
            push!(adj_w_i, p_ij / (1 - p_ij))
        end
        push!(adj_w_2, adj_w_i)
    end
    return adj_w_2
end

# Reweight the X decoding graph for correlated MWPM (hard-decision version).
# Uses the MWPM Z correction as a binary estimate of Z errors rather than
# the soft marginal probabilities from the worm.
function reweight_XZ_MWPM(edge_prob_list, qubit_mapping_1, qubit_mapping_2, L, p, N_1, N_2)
    qubit_probs = zeros(L^2)
    for i in 1:N_1
        for j in 1:length(edge_prob_list[i])
            qubit_probs[qubit_mapping_1[i][j]] = edge_prob_list[i][j]
        end
    end

    adj_w_2 = []
    for i in 1:N_2
        adj_w_i = []
        for j in 1:length(qubit_mapping_2[i])
            p_q  = qubit_probs[qubit_mapping_2[i][j]]
            w_ij = 0.5*p_q + (1 - p_q)*(2p/3)/(1 - 2p/3)
            push!(adj_w_i, w_ij)
        end
        push!(adj_w_2, adj_w_i)
    end
    return adj_w_2
end

# =============================================================================
# Simulation Parameters
# =============================================================================

L_list  = [21]       #L must be **odd** for this code.
N_samp  = 10000      # worm MCMC samples per syndrome
N_real  = 1          # Monte Carlo realisations per (L, p)
t_auto  = 1000       # worm autocorrelation time (steps between samples)
p_list = [0.08] 


# =============================================================================
# Main Simulation Loop
# =============================================================================
N_parse = 5       # number of X↔Z parsing iterations
nx = pyimport("networkx")

for L in L_list

    # 1. Build X and Z decoding graphs from the RSC parity check matrices.
    # Adjacency matrix A = H·Hᵀ; A[i,j] counts shared qubits between checks
    # i and j, giving the multi-graph structure of the decoding graph.
    n = L^2
    x_checks, x_log_qubits = find_X_checks(L)
    z_checks                = find_Z_checks(L)
    H_x = find_check_mat(x_checks, n)
    H_z = find_check_mat(z_checks, n)
    adj_x = H_x * transpose(H_x)
    adj_z = H_z * transpose(H_z)
    adj_x_comp = adj_to_compact_adj(adj_x)
    adj_z_comp = adj_to_compact_adj(adj_z)

    qubit_mapping_x = edge_to_qubit(adj_x_comp, H_x, n)
    qubit_mapping_z = edge_to_qubit(adj_z_comp, H_z, n)
    qubit_edges_x   = qubits_to_edge(qubit_mapping_x, n)
    qubit_edges_z   = qubits_to_edge(qubit_mapping_z, n)

    N_nodes_x = length(adj_x_comp)
    N_nodes_z = length(adj_z_comp)
    X_log_ind = find_X_logicals(x_log_qubits, qubit_mapping_x)

    for p in p_list
        p_ = 2*p/3      # marginal single-qubit X (or Z) rate under depolarizing noise

        correl_MLD_list_parses = [[] for _ in 1:N_parse]
        MWPM_uncorrel_list = []
        MWPM_correl_list   = []
        uncorrel_MLD_list  = []

        for _ in 1:N_real

            # 2. Initialise uniform weights (independent noise assumption)
            adj_x_w = [fill(p_/(1-p_), length(adj_x_comp[i])) for i in 1:N_nodes_x]
            adj_z_w = [fill(p_/(1-p_), length(adj_z_comp[i])) for i in 1:N_nodes_z]

            # 3. Sample a depolarizing error and compute X and Z syndromes
            q_x, q_z = generate_error(p, n)
            X_error   = x_error_edges(q_x, adj_x_comp, qubit_edges_x)
            synd_x    = find_syndrome(qubit_mapping_x, q_x)
            synd_z    = find_syndrome(qubit_mapping_z, q_z)

            # 4. Uncorrelated MWPM baseline (X and Z decoded independently)
            Z_MWPM = find_MWPM(synd_z, adj_z_w, N_nodes_z, adj_z_comp)
            X_MWPM = find_MWPM(synd_x, adj_x_w, N_nodes_x, adj_x_comp)

            push!(MWPM_uncorrel_list,
                  logical_measure(symmetric_dif(X_error, X_MWPM), X_log_ind) == 1 ? 1 : 0)

            # 5. Correlated MWPM: reweight X graph using the Z MWPM hard decision.
            # Convert MWPM bonds to binary error estimates: σ = -1 → p = 0.5, σ = 1 → p = 0
            error_dist_MWPM_Z     = [((-Z_MWPM[i] .+ 1) / 2) for i in 1:N_nodes_z]
            adj_x_MWPM_reweighted = reweight_XZ_MWPM(error_dist_MWPM_Z, qubit_mapping_z,
                                                       qubit_mapping_x, L, p_, N_nodes_z, N_nodes_x)
            X_MWPM_correl = find_MWPM(synd_x, adj_x_MWPM_reweighted, N_nodes_x, adj_x_comp)

            push!(MWPM_correl_list,
                  logical_measure(symmetric_dif(X_error, X_MWPM_correl), X_log_ind) == 1 ? 1 : 0)

            # 6. Uncorrelated worm decoding baseline (uniform weights, no Z info).
            # The returned marginal distribution is discarded here.
            MLD_worm_result, _ = worm_decoder(adj_x_w, adj_x_comp, X_MWPM, N_samp, t_auto,
                                               qubit_mapping_x, qubit_edges_x, X_error, X_log_ind)
            push!(uncorrel_MLD_list, MLD_worm_result)

            # 7. Iterative X↔Z parsing: alternate between decoding Z and X with
            # the worm, passing soft marginal probabilities from each run to
            # reweight the other graph. Each iteration refines the correlated
            # estimate; adj_z_w is updated in place between parse steps.
            for i_parse in 1:N_parse
            
                # Decode Z with worm to obtain soft marginal Z error probabilities
                global Z_matching_loop      = [fill(1, length(adj_z_comp[i])) for i in 1:N_nodes_z]
                global error_distribution_Z = [fill(0, length(adj_z_comp[i])) for i in 1:N_nodes_z]

                for _ in 1:N_samp
                    global Z_matching_loop = generate_independent_sample(
                        Z_matching_loop, t_auto, adj_z_w, adj_z_comp,
                        deepcopy(Z_MWPM), qubit_mapping_z, qubit_edges_z
                    )
                    Z_matching_MLD = symmetric_dif(Z_MWPM, Z_matching_loop)
                    global error_distribution_Z = [
                        (error_distribution_Z[i] .+ (-Z_matching_MLD[i] .+ 1) / (2*N_samp))
                        for i in 1:N_nodes_z
                    ]
                end

                # Reweight X graph with soft Z marginals and decode X with worm
                adj_x_reweighted = reweight_XZ(error_distribution_Z, qubit_mapping_z,
                                                qubit_mapping_x, L, p_, N_nodes_z, N_nodes_x)
                correl_MLD_result, error_distribution_X = worm_decoder(
                    adj_x_reweighted, adj_x_comp, X_MWPM, N_samp, t_auto,
                    qubit_mapping_x, qubit_edges_x, X_error, X_log_ind
                )
                push!(correl_MLD_list_parses[i_parse], correl_MLD_result)

                # Reweight Z graph with soft X marginals for the next parse iteration
                adj_z_w = reweight_XZ(error_distribution_X, qubit_mapping_x,
                                       qubit_mapping_z, L, p_, N_nodes_x, N_nodes_z)
            end
        end

        p_logical_fail_data = DataFrame(
            MWPM_correlated_success            = MWPM_correl_list,
            MWPM_uncorrelated_success          = MWPM_uncorrel_list,
            Worm_uncorrelated_decoding_success = uncorrel_MLD_list,
            Worm_parse_1                       = correl_MLD_list_parses[1],
            Worm_parse_2                       = correl_MLD_list_parses[2],
            Worm_parse_3                       = correl_MLD_list_parses[3],
            Worm_parse_4                       = correl_MLD_list_parses[4],
            Worm_parse_5                       = correl_MLD_list_parses[5],
        )

    
    end
    
end
