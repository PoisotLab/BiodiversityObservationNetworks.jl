"""
    _build_kdtree(cpool::CandidatePool)

Build a `KDTree` from candidate coordinates for `O(n*log(n))` nearest-neighbor queries.
"""
_build_kdtree(cpool::CandidatePool) = KDTree(Float32.(cpool.coordinates))

"""
    _neighbor_order(tree, coords, i)

Return all neighbors of candidate `i` ordered by distance (excluding self).
"""
function _neighbor_order(tree::KDTree, coords::Matrix, i::Int)
    idxs, _ = knn(tree, coords[:, i], size(coords, 2), true)
    return idxs[2:end] # first is always self
end

"""
    _neighbor_map(tree, coords)

Return a dictionary with keys for all nodes `i` pointing to the list of indices sorted by distance to `i`.
"""
function _neighbor_map(tree::KDTree, coords::Matrix)
    return Dict([i=>_neighbor_order(tree, coords, i) for i in 1:size(coords,2)])
end





#=
_node_distance(a::Tuple, b::Tuple) = sqrt(sum(a .- b) .^ 2)
_node_distance(a::CartesianIndex, b::CartesianIndex) = sqrt(sum(Tuple(a - b) .^ 2))


getdistancematrix(rd::RasterDomain) = [_node_distance(i,j) for i in vec(CartesianIndices(rd.data)), j in vec(CartesianIndices(rd.data))]

getdistancematrix(bon::BiodiversityObservationNetwork) = [_node_distance(i,j) for i in bon.nodes, j in bon.nodes]


function getnearestneighbors(bon::BiodiversityObservationNetwork)
    coord_tuples = bon.nodes
    coord_mat = _get_coord_matrix(vec([x for x in bon]))
    return _create_nearest_neighbor_map(coord_tuples, coord_mat)
end

function getnearestneighbors(rd::RasterDomain)
    indices = getpool(rd)
    flat_indices = Vector(vec(indices))
    coords = _get_coord_matrix(flat_indices)
    return _create_nearest_neighbor_map(flat_indices, coords)
end

function _create_nearest_neighbor_map(coord_vec::Vector{T}, coord_mat) where {T}
    tree = KDTree(coord_mat)
    k = length(coord_vec)
    neighbor_map = Dict{T,Vector{T}}()
    idx_map = Dict()
    for (idx, ci) in enumerate(coord_vec)
        neighbor_ids, _ = knn(tree, coord_mat[:, idx], k, true)

        # first idx will be self 
        neighbor_map[ci] = coord_vec[neighbor_ids][2:end]
        idx_map[idx] = neighbor_ids[2:end]
    end
    return neighbor_map, idx_map
end

function _get_coord_matrix(vec::Vector{<:Tuple})
    F = Sys.WORD_SIZE == 64 ? Float64 : Float32
    coord_mat = F.(hcat([[x...] for x in vec]...))  # each column is a point
    return coord_mat
end

function _get_coord_matrix(vec::Vector{<:CartesianIndex})
    F = Sys.WORD_SIZE == 64 ? Float64 : Float32
    coord_mat = F.(hcat([[ci...] for ci in Tuple.(vec)]...))  # each column is a point
    return coord_mat
end

=#