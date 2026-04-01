"""
    CandidatePool{K}

The internal representation of a sampling domain: a matrix of `n` candidate locations, each with coordinates. Each candidate coordinate can potentially have features (auxiliary variables) and inclusion weights, which may be required/optional for some sampling methods. 

# Fields
- `n::Int` — the number of candidate locations
- `coordinates::Matrix` — `2 x n` spatial coordinate matrix
- `features::Union{Matrix, Missing}` — `p x n` feature matrix (where p is the number of auxiliary variables) or `Missing`
- `inclusion::Vector` — `n`-vector of relative inclusion probability summing to 1
- `keys::Vector{K}` — original identifiers (e.g. `CartesianIndex` for raster inputs)
- `metadata::Dict{Symbol, Any}` — source info, CRS, extent, etc.
"""
struct CandidatePool{K}
    n::Int
    coordinates::Matrix
    features::Union{Matrix, Missing}
    inclusion::Vector
    keys::Vector{K}
    metadata::Dict{Symbol, Any}
end

# Base overloads 
Base.length(cp::CandidatePool) = cp.n
Base.size(cp::CandidatePool) = (cp.n,)
Base.iterate(cp::CandidatePool, state = 1) = iterate(cp.keys, state)
Base.getindex(cp::CandidatePool, i::Int) = cp.keys[i]

# ========================================================================
# Constructors
# ========================================================================
