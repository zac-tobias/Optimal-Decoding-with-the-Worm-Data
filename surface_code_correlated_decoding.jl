using Random
using LinearAlgebra
using CSV
using DataFrames
using Graphs
using Random
using PyCall

function find_X_checks(L)
    x_checks_list = []
    x_logical  = []
    for i in 1:(L-1)/2
        x = Int(2*i-1)

        for j in 1:(L-1)/2
            y = Int(2*j-1)

            x_check_1 = [coords_to_qubit(x+1,y,L),coords_to_qubit(x+2,y,L),coords_to_qubit(x+1,y+1,L),coords_to_qubit(x+2,y+1,L)]
            x_check_2 = [coords_to_qubit(x,y+1,L),coords_to_qubit(x,y+2,L),coords_to_qubit(x+1,y+1,L),coords_to_qubit(x+1,y+2,L)]

            push!(x_checks_list,x_check_1)
            push!(x_checks_list,x_check_2)
        end
    end
    for k in 1:(L-1)/2
        x = Int(2*k-1)
        x_boundary_check_1 = [coords_to_qubit(x,1,L),coords_to_qubit(x+1,1,L)]
        x_boundary_check_2 = [coords_to_qubit(x+1,L,L),coords_to_qubit(x+2,L,L)]
        push!(x_checks_list, x_boundary_check_1)
        push!(x_checks_list, x_boundary_check_2)
    end
    boundary_node_x = []
    for y in 1:L
        push!(boundary_node_x,coords_to_qubit(1,y,L))
        push!(x_logical,coords_to_qubit(1,y,L))
        push!(boundary_node_x,coords_to_qubit(L,y,L))
    end
    push!(x_checks_list,boundary_node_x)
    return(x_checks_list,x_logical)

end

function find_Z_checks(L)
    z_checks_list = []
    for i in 1:(L-1)/2
        x = Int(2*i-1)

        for j in 1:(L-1)/2
            y = Int(2*j-1)

            z_check_1 = [coords_to_qubit(x,y,L),coords_to_qubit(x+1,y,L),coords_to_qubit(x,y+1,L),coords_to_qubit(x+1,y+1,L)]
            z_check_2 = [coords_to_qubit(x+1,y+1,L),coords_to_qubit(x+1,y+2,L),coords_to_qubit(x+2,y+1,L),coords_to_qubit(x+2,y+2,L)]

            push!(z_checks_list,z_check_1)
            push!(z_checks_list,z_check_2)
        end
    end
    for k in 1:(L-1)/2
        y = Int(2*k-1)
        z_boundary_check_1 = [coords_to_qubit(L,y,L),coords_to_qubit(L,y+1,L)]
        z_boundary_check_2 = [coords_to_qubit(1,y+1,L),coords_to_qubit(1,y+2,L)]
        push!(z_checks_list,  z_boundary_check_1)
        push!(z_checks_list, z_boundary_check_2)
    end

    boundary_node_z = []
    for x in 1:L
        push!(boundary_node_z,coords_to_qubit(x,1,L))
        push!(boundary_node_z,coords_to_qubit(x,L,L))
    end
    push!(z_checks_list,boundary_node_z)

    return(z_checks_list)
end

function coords_to_qubit(x,y,L)
    return(x+(y-1)*L)
end

function find_check_mat(checks_list,n)
    m = length(checks_list)
    H = zeros(Int, m, n)
    #println(m)
    for i in 1:m 
        check_m = checks_list[i]
        #println(check_m)
        for qu in check_m
            #println(qu)
            H[i,qu] = 1
        end
    end
    return(H)
end

function adj_to_compact_adj(adj)
    N = length(adj[1,:])
    adj_compact = []
    for i in 1:N 
        adj_i = adj[i,:]
        adj_i_compact = []
        for k in 1:N 
            adj_ik = adj_i[k]
            if adj_ik > 0 
                if k != i
                    for _ in 1:adj_ik
                        push!(adj_i_compact,k)
                    end
                end
            end
        end
        push!(adj_compact,adj_i_compact)
    end
    return(adj_compact)
end

function edge_to_qubit(adj_comp,H,n)
    qubit_list_mapping = []
    N = length(adj_comp)
    H_q = transpose(H)
    for m_1 in 1:N
        adj_comp_i = adj_comp[m_1]
        qubit_list_i = []
        for m_2 in adj_comp_i
            for i in 1:n
                
                h_q_i = H_q[i,:]
                h_q_i_ind = findall(x -> x == 1, h_q_i)
                
                if m_1 in h_q_i_ind
                    if m_2 in h_q_i_ind
                        if i ∉ qubit_list_i
                            q = i
                            push!(qubit_list_i, q)
                        end
                    end
                end
            end
        
        end
        push!(qubit_list_mapping,qubit_list_i)
    end
    return(qubit_list_mapping)
end

function edge_probs_to_qubits(edge_prob_list,qubit_mapping)
    qubit_probs = [0.0 for _ in 1:n]
    println(qubit_probs)
    N = length(edge_prob_list)    
    for i in 1:N
        N_i = length(edge_prob_list[i])
        for j in 1:N_i
            p_ij = edge_prob_list[i][j]
            q_ij = qubit_mapping[i][j]
            #println(q_ij)
            qubit_probs[q_ij] = p_ij
        end
    end
    return(qubit_probs)
end

function qubits_to_edge_prob(q_prob_list,qubit_mapping)
    edge_probs = []
    for edges_i in qubit_mapping
        p_list_i = []
        for qubit_ij in edges_i
            push!(p_list_i,q_prob_list[qubit_ij])
        end
        push!(edge_probs,p_list_i)
    end
    return(edge_probs)
end

function generate_error(p, n)
    qubit_list_x = Int[]
    qubit_list_z = Int[]

    for i in 1:n
        x = rand()
        if x < p/3
            push!(qubit_list_x, 1)
            push!(qubit_list_z, 1)
        elseif x < 2p/3
            push!(qubit_list_x, 1)
            push!(qubit_list_z, 0)
        elseif x < p
            push!(qubit_list_x, 0)
            push!(qubit_list_z, 1)
        else
            push!(qubit_list_x, 0)
            push!(qubit_list_z, 0)
        end
    end

    return qubit_list_x, qubit_list_z
end

function find_syndrome(qubit_mapping,qubits)
    m = length(qubit_mapping)
    syndrome = []
    
    for i in 1:m
        qubit_mapping_i = qubit_mapping[i]
        parity = 1
        for q_ind in qubit_mapping_i
            
            if qubits[q_ind] == 1
              
                parity = -parity
            end
        end
        if parity == -1
            push!(syndrome,i)
        end
    end
    return(syndrome)
end

function generate_independent_sample(original_loop_config,t_auto,A_weights, Adj_matrix, Quenched_Disorder,qubit_mapping,qubit_edges)
    loop = deepcopy(original_loop_config)
    c = 0
    while c < t_auto
        new_loop = generate_new_loop_config(loop,A_weights,Adj_matrix,Quenched_Disorder,qubit_mapping,qubit_edges)
        loop = new_loop
        c += 1
    end
    return(loop)
end

function generate_new_loop_config(loop_config, Weights_Matrix, Adjacency_matrix, quenched_dis_matrix,edge_to_qubit,qubit_to_edge)
    new_loop_config = false
    N_sites = length(Adjacency_matrix)
    start_site = rand(1:N_sites)
    current_site = copy(start_site)
    while !new_loop_config
        direction = rand(1:length(Adjacency_matrix[current_site]))
        new_site = Adjacency_matrix[current_site][direction]
        J_ij = quenched_dis_matrix[current_site][direction]
        sigma_ij = loop_config[current_site][direction]
        if Int(0.5 *(1+J_ij * sigma_ij)) == 0
            w_ij = (1/(Weights_Matrix[current_site][direction])) * (length(Weights_Matrix[current_site])/length(Weights_Matrix[new_site]))
        else
            w_ij = Weights_Matrix[current_site][direction] * (length(Weights_Matrix[current_site])/length(Weights_Matrix[new_site]))
        end

        if w_ij > 1
            w_ij = 1
        end

        r = rand()
        if r < w_ij
            loop_config[current_site][direction] = -sigma_ij
            q_ij = edge_to_qubit[current_site][direction]

            qubit_edges_ij = qubit_to_edge[q_ij]
            qubit_edges_1 = qubit_edges_ij[1]
            qubit_edges_2 = qubit_edges_ij[2]
            #println(new_site)
            #println(qubit_edges_1[1])
            #println(qubit_edges_2[1])
            #println("-------")
            if qubit_edges_1[1] == new_site
                old_index = qubit_edges_1[2]
            else
                old_index = qubit_edges_2[2]
            end
            
            loop_config[new_site][old_index] = -sigma_ij
            current_site = new_site

            #qubit = edge_to_qubit[current_site][direction]
            

        end
        if current_site == start_site
                new_loop_config = true
        end

    end

    return loop_config
end

function find_MWPM(synd,Adj_w,N_nodes,adj)
    Gr = nx.Graph()

    for (i, nbrs) in enumerate(adj)
        for (k, j) in enumerate(nbrs)     # k is index into adj_w[i]
            if j > i                      # avoid double-adding edges
                w = -log(Adj_w[i][k])           # get the weight from adj_w
                Gr.add_edge(i, j, weight=w)
            end
        end
    end

    CG = nx.Graph()
    for (idx, u) in enumerate(synd)
        for v in synd[idx+1:end]
            dist = nx.shortest_path_length(Gr, u, v, weight="weight")
            CG.add_edge(u, v, weight = -dist)
        end
    end

    # --- Run Blossom MWPM on complete graph ---
    matching = nx.max_weight_matching(CG, maxcardinality=true)
    # --- Recover the actual paths in the original graph ---
    
    path_list = []
    for (u,v) in matching
        path = nx.shortest_path(Gr, u, v, weight="weight")
        for i in 1:length(path)-1
            push!(path_list, [path[i], path[i+1]])
        end
        edges = [[path[i], path[i+1]] for i in 1:length(path)-1]
        
    end

    matched_config =  [fill(1, length(adj[i])) for i in 1:N_nodes]

    for edge_ij in path_list
        node_i = edge_ij[1]
        node_j = edge_ij[2]
        ind_ij = findfirst( x -> x == node_j, adj[node_i])
        ind_ji = findfirst( x-> x == node_i, adj[node_j])

        matched_config[node_i][ind_ij] = -1
        matched_config[node_j][ind_ji] = -1

    end

    return(matched_config)
end

function find_matching_any(synd_list, Adj_w, N_nodes, Adj)
    Gr = nx.Graph()

    for (i, nbrs) in enumerate(Adj)
        for (k, j) in enumerate(nbrs)
            if j > i && Adj_w[i][k] > 0
                w = -log(Adj_w[i][k])
                Gr.add_edge(i, j, weight=w)
            end
        end
    end

    n_synd = Int(length(synd_list) ÷ 2)
    matched_config_full = [fill(1, length(Adj[i])) for i in 1:N_nodes]

    for i in 1:n_synd
        CG = nx.Graph()
        synd = [synd_list[2*i-1], synd_list[2*i]]
        u, v = synd
        dist = nx.shortest_path_length(Gr, u, v, weight="weight")
        CG.add_edge(u, v, weight=-dist)

        matching = nx.max_weight_matching(CG, maxcardinality=true)

        path_list = []
        for (u, v) in collect(matching)
            path = nx.shortest_path(Gr, u, v, weight="weight")
            for j in 1:length(path)-1
                push!(path_list, [path[j], path[j+1]])
            end
        end

        matched_config = [fill(1, length(Adj[i])) for i in 1:N_nodes]
        for edge_ij in path_list
            node_i, node_j = edge_ij
            ind_ij = findfirst(x -> x == node_j, Adj[node_i])
            ind_ji = findfirst(x -> x == node_i, Adj[node_j])
            matched_config[node_i][ind_ij] = -1
            matched_config[node_j][ind_ji] = -1
        end

        matched_config_full = symmetric_dif(matched_config_full, matched_config)
    end

    return matched_config_full
end

function symmetric_dif(a1,a2)
    a_12 = []
    for j in 1:length(a1)
        a_12_i = []
        for i in 1:length(a1[j])
            a_12_ij = a1[j][i] * a2[j][i]
            push!(a_12_i,a_12_ij)
        end
        push!(a_12,a_12_i)
    end
    return(a_12)
end

function find_syndrome_edge_conf(dis)
    synd_list = []
    for i in 1:length(dis)
        dis_i = dis[i]
        l = 1
        for dis_ij in dis_i
            l = l*dis_ij
        end
        
        if l == -1
            push!(synd_list,i)
        end
    end
    return(synd_list)
end

function qubits_to_edge(qubit_mapping, n_qubits)
    qubit_edges = [ [] for i in 1:n_qubits]
    n_edges = length(qubit_mapping)
    for i in 1:n_edges
        qubit_mapping_i = qubit_mapping[i]
        for j in 1:length(qubit_mapping_i)
            q_ij = qubit_mapping_i[j]
            push!(qubit_edges[q_ij],[i,j])
        end
    end
    return(qubit_edges)
end

function reweight_XZ(edge_prob_list, qubit_mapping_1,qubit_mapping_2,L,p,N_x,N_z)
    qubit_probs = zeros(L^2)
    for i in 1:N_x
        edge_prob_list_i = edge_prob_list[i]
        for j in 1:length(edge_prob_list_i)
            q_ij = qubit_mapping_1[i][j]
            qubit_probs[q_ij] = edge_prob_list_i[j]
        end
    end
    
    adj_w_2 = []


    for i_ in 1:N_z
        adj_w_i = []
        q_mapping_i = qubit_mapping_2[i_]
        for j in 1: length(q_mapping_i)
            q_ij = q_mapping_i[j]
            p_q_ij = qubit_probs[q_ij]

            p_ij = 0.5*p_q_ij + (1-p_q_ij)*(p/3)/(1-2*p/3)
            
            push!(adj_w_i, p_ij/(1-p_ij))
        end
        push!(adj_w_2,adj_w_i)
    end
     
    #println(adj_w_2)
    return(adj_w_2) 
end
    
function find_X_logicals(x_log_qubit, edge_to_qubits)
    x_log_ind = []
    N_x = length(edge_to_qubits)
    boundary_qubits = edge_to_qubits[N_x]
    for q in x_log_qubit
        if q in boundary_qubits
            q_ind  = findfirst(==(q), boundary_qubits) 
            push!(x_log_ind,q_ind)
        end
    end
    return(x_log_ind)
end

function x_error_edges(q_x_errors, adj_x,q_edge_map)
    N_x = length(adj_x)
    x_error_edge_config =[fill(1, length(adj_x[i])) for i in 1:N_x]
    #println(q_x_errors)
    
    for i in 1:length(q_x_errors)
        q = q_x_errors[i]
        if q ==1
            q_edges = q_edge_map[i]
            q_edges_1 = q_edges[1]
            q_edges_2 = q_edges[2]
            #println(q_edges_1)
            x_error_edge_config[q_edges_1[1]][q_edges_1[2]] = -1
            x_error_edge_config[q_edges_2[1]][q_edges_2[2]] = -1
        end
    end

    return(x_error_edge_config)
end

function logical_measure(X_config, log_indices)
    l = 1
    N_x_ = length(X_config)
    X_conf_boundary = X_config[N_x_]
    for i in log_indices
        l = l* X_conf_boundary[i]
    end
    return(l)
end

function worm_decoder(adj_w,adj,matching,n_samples,t_auto,edges_to_qubits,quibit_edges,error_chain,log_ind)
    sector_tally = [0,0]
    matching_loop = [fill(1, length(adj[i])) for i in 1:length(adj)]
    for _ in 1:n_samples
        matching_loop = generate_independent_sample(matching_loop,t_auto, adj_w, adj,deepcopy(matching),edges_to_qubits,quibit_edges)
        matching_MLD = symmetric_dif(matching,matching_loop)

        Total_op = symmetric_dif(error_chain,matching_MLD)
        
        logical_val = logical_measure(Total_op,log_ind)
        
        if logical_val == 1
            sector_tally[1]+=1
        else
            sector_tally[2]+=1
        end
    end

    if sector_tally[1] == sector_tally[2]
        r = rand([0, 1])
        return(r)
    elseif sector_tally[1] > sector_tally[2]
        return(1)
    else
        return(0)
    end

end

function reweight_XZ(edge_prob_list, qubit_mapping_1,qubit_mapping_2,L,p,N_x,N_z)
    qubit_probs = zeros(L^2)
    for i in 1:N_x
        edge_prob_list_i = edge_prob_list[i]
        for j in 1:length(edge_prob_list_i)
            q_ij = qubit_mapping_1[i][j]
            qubit_probs[q_ij] = edge_prob_list_i[j]
        end
    end
    
    adj_w_2 = []


    for i_ in 1:N_z
        adj_w_i = []
        q_mapping_i = qubit_mapping_2[i_]
        for j in 1: length(q_mapping_i)
            q_ij = q_mapping_i[j]
            p_q_ij = qubit_probs[q_ij]

            p_ij = 0.5*p_q_ij + (1-p_q_ij)*(p/3)/(1-2*p/3)
            
            push!(adj_w_i, p_ij/(1-p_ij))
        end
        push!(adj_w_2,adj_w_i)
    end
     
    #println(adj_w_2)
    return(adj_w_2) 
end

function reweight_XZ_MWPM(edge_prob_list, qubit_mapping_1,qubit_mapping_2,L,p,N_x,N_z)
    qubit_probs = zeros(L^2)
    for i in 1:N_x
        edge_prob_list_i = edge_prob_list[i]
        for j in 1:length(edge_prob_list_i)
            q_ij = qubit_mapping_1[i][j]
            qubit_probs[q_ij] = edge_prob_list_i[j]
        end
    end
    
    adj_w_2 = []


    for i_ in 1:N_z
        adj_w_i = []
        q_mapping_i = qubit_mapping_2[i_]
        for j in 1: length(q_mapping_i)
            q_ij = q_mapping_i[j]
            p_q_ij = qubit_probs[q_ij]
            
            w_ij = 0.5*p_q_ij + (1-p_q_ij)*(2*p/3)/(1-2*p/3)
            
            push!(adj_w_i,w_ij)
        end
        push!(adj_w_2,adj_w_i)
    end
     
    #println(adj_w_2)
    return(adj_w_2) 
end


println("--------------------------------------------------")
L_list = [27]
N_samp = 2000
N_real = 50
#t_auto = 100
t_auto = 300
#p_list = [0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2]
#p_list = [0.17,0.18,0.19,0.2]
#p_list = [0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2]
p_list = [0.17,0.19]
#---------------------------------
nx = pyimport("networkx")
for L in L_list
    
    n = L^2
    x_checks, x_log_qubits = find_X_checks(L)
    z_checks = find_Z_checks(L)
    H_x = find_check_mat(x_checks,n)
    H_z = find_check_mat(z_checks,n)
    adj_x = H_x * transpose(H_x)
    adj_z =  H_z * transpose(H_z)
    adj_x_comp = adj_to_compact_adj(adj_x)
    adj_z_comp = adj_to_compact_adj(adj_z)
    qubit_mapping_x = edge_to_qubit(adj_x_comp,H_x,n)
    qubit_mapping_z = edge_to_qubit(adj_z_comp,H_z,n)
    qubit_edges_x  = qubits_to_edge(qubit_mapping_x,n)
    qubit_edges_z  = qubits_to_edge(qubit_mapping_z,n)
    N_nodes_x = length(adj_x_comp)
    N_nodes_z = length(adj_z_comp)
    X_log_ind = find_X_logicals(x_log_qubits, qubit_mapping_x)
     
    
     
    for p in p_list
        println(p)
        p_ = 2*p/3
    
        correl_MLD_list = []
        MWPM_uncorrel_list = []
        MWPM_correl_list = []
        uncorrel_MLD_list = []
    
        adj_x_w = [fill(p_/(1-p_), length(adj_x_comp[i])) for i in 1:N_nodes_x]
        adj_z_w = [fill(p_/(1-p_), length(adj_z_comp[i])) for i in 1:N_nodes_z]
    
        for _ in 1:N_real
    
            q_x, q_z = generate_error(p,n)
            X_error = x_error_edges(q_x,adj_x_comp,qubit_edges_x)
    
    
            synd_x = find_syndrome(qubit_mapping_x,q_x)
            synd_z = find_syndrome(qubit_mapping_z,q_z)
        
    
            Z_MWPM = find_MWPM(synd_z,adj_z_w,N_nodes_z,adj_z_comp)
            X_MWPM = find_MWPM(synd_x,adj_x_w,N_nodes_x,adj_x_comp)
    
            MWPM_total_op = symmetric_dif(X_error,X_MWPM)
            x_log_MWPM = logical_measure(MWPM_total_op,X_log_ind)
            if x_log_MWPM == 1
                push!(MWPM_uncorrel_list,1)
            else
                push!(MWPM_uncorrel_list,0)
            end
            error_dist_MWPM_Z =  [((-Z_MWPM[i].+1)/(2)) for i in 1:N_nodes_z]
            adj_x_MWPM_reweighted = reweight_XZ_MWPM(error_dist_MWPM_Z,qubit_mapping_z,qubit_mapping_x,L,p_,N_nodes_z,N_nodes_x)
            X_MWPM_correl = find_MWPM(synd_x,adj_x_MWPM_reweighted,N_nodes_x,adj_x_comp)
    
            MWPM_correl_total_op = symmetric_dif(X_error,X_MWPM_correl)
            x_log_MWPM_correl = logical_measure(MWPM_correl_total_op,X_log_ind)
            if x_log_MWPM_correl == 1
                push!(MWPM_correl_list,1)
            else
                push!(MWPM_correl_list,0)
            end
    
    
    
            global Z_matching_loop = [fill(1, length(adj_z_comp[i])) for i in 1:N_nodes_z]
            global error_distribution_Z = [fill(0, length(adj_x_comp[i])) for i in 1:N_nodes_z]
            
    
            for _ in 1:N_samp
                global Z_matching_loop = generate_independent_sample(Z_matching_loop,t_auto,adj_z_w,adj_z_comp,deepcopy(Z_MWPM),qubit_mapping_z,qubit_edges_z)
                #println(find_syndrome_edge_conf(X_matching_loop))
                #println(adj_w)
                Z_matching_MLD = symmetric_dif(Z_MWPM,Z_matching_loop)
                global error_distribution_Z = [(error_distribution_Z[i] .+ (-Z_matching_MLD[i].+1)/(2*N_samp)) for i in 1:N_nodes_z]
            end
    
            #println("---------------------------")
            adj_x_reweighted = reweight_XZ(error_distribution_Z,qubit_mapping_z,qubit_mapping_x,L,p_,N_nodes_z,N_nodes_x)
            correl_MLD_result = worm_decoder(adj_x_reweighted,adj_x_comp,X_MWPM,N_samp,t_auto,qubit_mapping_x,qubit_edges_x,X_error,X_log_ind)
            push!(correl_MLD_list,correl_MLD_result)
            MLD_worm_result = worm_decoder(adj_x_w,adj_x_comp,X_MWPM,N_samp,t_auto,qubit_mapping_x,qubit_edges_x,X_error,X_log_ind)
            push!(uncorrel_MLD_list,MLD_worm_result)
    
        end
    

    folder_path = "/home/users/ztobias/Decoding_Project_Data/Data_for_paper/Correlated_worm_decoding"
    filename = joinpath(
                folder_path,
                "Correlated_worm_decoding_RSC_L_" * string(L) * "_p_" * string(round(Int, 10000 * p)) * ".csv"
            )
            
    col_names = ["MWPM_correlated_success", "MWPM_uncorrelated_success", "Worm_correlated_decoding_success","Worm_uncorrelated_decoding_success"]
            
            # Create file with header if it doesn't exist
    if !isfile(filename)
            open(filename, "w") do io
                write(io, join(col_names, ",") * "\n")
        end
    end
            
            # Construct the DataFrame from the lists
    p_logical_fail_data = DataFrame(
                MWPM_correlated_success = MWPM_correl_list,
                MWPM_uncorrelated_success = MWPM_uncorrel_list,
                Worm_correlated_decoding_success = correl_MLD_list,
                Worm_uncorrelated_decoding_success = uncorrel_MLD_list
            )
            
            # Append the data to the file
            
    CSV.write(filename, p_logical_fail_data; append=true, writeheader=false)
       
    
    end
end
