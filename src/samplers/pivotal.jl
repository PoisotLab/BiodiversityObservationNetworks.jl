Base.@kwdef struct Pivotal{I<:Integer} <: BONSampler
    number_of_nodes::I = 100
    maximum_iterations::I = 10^6
end 

_valid_geometries(::Pivotal) = (BiodiversityObservationNetwork)


# 1. Randomly choose a unit i, and a nearest neighbor to i, j. (If two or more
# units are equidistant, choose at random)
# 2. Update inclusion πᵢ and πⱼ with the following rule:
#       
#  If πᵢ + πⱼ < 1:
#       (πᵢ', πⱼ') = (0, πᵢ + πⱼ)   w.p. πⱼ / (πᵢ + πⱼ) and 
#       (πᵢ + πⱼ, 0)   w.p. πᵢ / (πᵢ + πⱼ)  -- (same as 1 - above prob)
#  Else:   
#       (πᵢ', πⱼ') = (1, πᵢ + πⱼ - 1)   w.p. (1 - πⱼ) / (2 - πᵢ - πⱼ)  
#       (πᵢ', πⱼ') = (πᵢ + πⱼ - 1, 1)   w.p. (1 - πᵢ) / (2 - πᵢ - πⱼ) -- (same as 1 - above prob) 
# 
# 3. Repeat (1) and (2) until all inclusion probabilities are 0 or 1
# 

function _get_distance_matrix(bon::BiodiversityObservationNetwork)
    n_nodes = size(bon)
    dist_mat = zeros(n_nodes, n_nodes)
    
    for i in 1:n_nodes, j in 1:n_nodes
        dist_mat[i,j] = i == j ? Inf : sqrt(sum((bon[i].coordinate .- bon[i].coordinate).^2))        
    end
    return dist_mat
end

function _pick_nearest(distance_matrix, idx)
    _, i = findmin(distance_matrix[idx,:])
    return i
end

function _mark_as_complete!(complete_idx, distance_matrix, idx)
    complete_idx[idx] = 1
    distance_matrix[idx, :] .= Inf
    distance_matrix[:, idx] .= Inf
end

function _above_one_update!(Π, complete_idx, distance_matrix, i, j)
    πᵢ, πⱼ = Π[i], Π[j] 
    add_idx, sub_idx = rand() < (1 - πⱼ)/(2 - πᵢ - πⱼ) ? (i, j) : (j, i)

    Π[add_idx] = 1     
    Π[sub_idx] = πᵢ + πⱼ - 1
    _mark_as_complete!(complete_idx, distance_matrix, add_idx)
end

function _below_one_update!(Π, complete_idx, distance_matrix, i, j)
    πᵢ, πⱼ = Π[i], Π[j] 
    add_idx, sub_idx = rand() < (πⱼ / (πᵢ + πⱼ)) ? (i, j) : (j, i)

    Π[add_idx] = πᵢ + πⱼ 
    Π[sub_idx] = 0
    _mark_as_complete!(complete_idx, distance_matrix, sub_idx)
end 

function _apply_update_rule!(Π, complete_idx, distance_matrix, i, j)
    if Π[i] + Π[j]  < 1 
        _below_one_update!(Π, complete_idx, distance_matrix, i, j)
    else
        _above_one_update!(Π, complete_idx, distance_matrix, i, j)
    end 
end 

function _sample(sampler::Pivotal, bon::BiodiversityObservationNetwork) 
    distance_matrix = _get_distance_matrix(bon)
    complete_idx = zeros(Bool, size(bon))

    base_π = sampler.number_of_nodes / size(bon)
    Π = zeros(size(bon))
    Π .= base_π

    ct = 0
    while !all(complete_idx)
        candidate_i = findall(!isone, complete_idx)
        i = rand(candidate_i)
        j = _pick_nearest(distance_matrix, i)

        _apply_update_rule!(Π, complete_idx, distance_matrix, i, j)

        ct += 1
        ct > sampler.maximum_iterations && break
    end 
    BiodiversityObservationNetwork(bon[findall(isone, Π)])
end 


