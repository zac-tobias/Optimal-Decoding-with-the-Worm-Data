using JSON
using Random
using LinearAlgebra
using DataFrames
using Graphs
using PyCall


function generate_new_loop_config(loop_config, Weights_Matrix, Adjacency_matrix, quenched_dis_matrix)
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
            old_index = findfirst(==(current_site), Adjacency_matrix[new_site])            
            loop_config[new_site][old_index] = -sigma_ij
            current_site = new_site
        end
        if current_site == start_site
            new_loop_config = true
        end
        
    end
    
    return loop_config
end

function topological_sector(loop_config, logical_pairs_list,adj)
    #Note this needs to be generalized for more than one logical
    k = 1
    logical_state = []
    for i in 1:k
        l_i_status = 1
        logical_i = logical_pairs_list
        for l_i_edge in logical_i
            l_i_ind = findfirst(==(l_i_edge[2]),adj[l_i_edge[1]])
            l_i_status = l_i_status*loop_config[l_i_edge[1]][l_i_ind]            
        end
        push!(logical_state,l_i_status)
    end

    return logical_state
end

function generate_independent_sample(original_loop_config,t_auto,A_weights, Adj_matrix, Quenched_Disorder)
    loop = copy(original_loop_config)
    c = 0
    while c < t_auto
        new_loop = generate_new_loop_config(loop,A_weights,Adj_matrix,Quenched_Disorder)
        loop = new_loop
        c += 1
    end
    return(loop)
end    

function run_sampling(initial_loop, N_samples, t_THERMALIZE, t_AUTO, A_w, A_, Quenched_Disorder, logical_pairs_list)
    loop = generate_independent_sample(initial_loop, t_THERMALIZE, A_w,  A_, Quenched_Disorder)
    logical_state_counter = [0,0]
    for _ in 1:N_samples
        new_loop = generate_independent_sample(loop, t_AUTO, A_w,  A_, Quenched_Disorder)
        logical_states = topological_sector(new_loop, logical_pairs_list,A_)
        loop = new_loop
        #Note this needs to be generalized for when there are multiple logicals
        if logical_states[1] == 1
            logical_state_counter[1] +=1
        else 
            logical_state_counter[2] += 1
        end        
       
    end
    return logical_state_counter
end


function Adj_compressed_to_A_mat(Adj_c,N_nodes)
    A_mat = zeros(Int, N_nodes,N_nodes)
    for i in 1:N_nodes
        Adj_i = copy(Adj_c[i])   
        for k in 1:length(Adj_i)
            A_mat[i,Adj_i[k]] = 1
        end
    end
    return(A_mat)
end

function find_MWPM(synd,Adj_w,N_nodes,Adj)
    Gr = nx.Graph()
    for (i, nbrs) in enumerate(Adj)
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

    matching = nx.max_weight_matching(CG, maxcardinality=true)
    
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



#---------------------------
#Inputs from STIM
filename_DEM = ARGS[1]
filename_synd = ARGS[2]

# --- Load data ---
data_DEM = JSON.parsefile(filename_DEM)
data_synd = JSON.parsefile(filename_synd)

# --- Extract fields ---
logicals       = data_DEM["logicals"]
adj            = data_DEM["adj"]
adj_w          = data_DEM["adj_w"]
t_AUTO         = data_DEM["t_autocorrelation"]
t_THERMALIZE   = data_DEM["t_thermalize"]
N_samples      = data_DEM["N_samples"]

syndrome_list  = data_synd["syndromes"]
rm(filename_synd; force=true)
#-------------------------------------
 
#Outputs
#Logical change from error

#Programme starts here
#-----------------------------
N_nodes = length(adj)
nx = pyimport("networkx")
correction_output_list = []
correction_prob_list = []
 

for syndrome in syndrome_list  
    trial_matching_ =  find_MWPM(syndrome,adj_w,N_nodes,adj)
    initial_log = topological_sector(trial_matching_,logicals,adj)
    initial_loop = [fill(1, length(adj[i])) for i in 1:N_nodes] 
    logical_sampling_results = run_sampling(initial_loop,N_samples,t_THERMALIZE,t_AUTO,adj_w,adj,trial_matching_,logicals)    

    if logical_sampling_results[1] > N_samples/2
        l = 1
        push!(correction_prob_list, logical_sampling_results[1]/N_samples)
    else
        l = -1
        push!(correction_prob_list, logical_sampling_results[2]/N_samples)
    end
       
    push!(correction_output_list,l * initial_log[1])
end




