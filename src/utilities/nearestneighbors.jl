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



