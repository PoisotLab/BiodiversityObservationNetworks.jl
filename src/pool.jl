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
# Helpers
# ========================================================================

function _get_coordinates_from_cartesian_indices(indices::Vector{CartesianIndex{2}})
    return hcat([[x[1], x[2]] for x in indices]...)
end 

_process_inclusion(::Missing, n_candidates::Int) = fill(1.0 / n_candidates, n_candidates)
function _process_inclusion(weights::AbstractVector, n_candidates::Int)
    length(weights) == n_candidates || throw(
        ArgumentError("Inclusion vector length ($(length(raw))) must match candidate count ($n_candidates)"))
    all(weights .>= 0) || throw(ArgumentError("Inclusion weights must be non-negative"))

    s = sum(weights)
    s > 0 || throw(ArgumentError("Inclusion weights must have positive sum"))

    weights ./= s
    return weights
end

_extract_and_process_inclusion(::Missing, keys, n) = _process_inclusion(missing, n)
function _extract_and_process_inclusion(inclusion::AbstractMatrix, keys, n)
    _process_inclusion([inclusion[k] for k in keys], n)
end


# ========================================================================
# Constructors
# ========================================================================

"""
    CandidatePool(mat::AbstractMatrix; mask = missing, inclusion = missing)

Convert a matrix into a [`CandidatePool`](@ref). Each cell becomes a candidate
with grid-index coordinates `(row, col)`. An optional boolean `mask` restricts
which cells are included. An optional `inclusion` matrix provides per-cell weights.
"""
function CandidatePool(mat::AbstractMatrix; mask = missing, inclusion = missing)
    # Check input validity
    ismissing(mask) || size(mask) == size(mat) ||
        throw(ArgumentError("Mask size $(size(mask)) must match domain size $(size(mat))"))

    # Apply mask
    all_indices = CartesianIndices(mat)
    valid = ismissing(mask) ? trues(size(mat)) : BitMatrix(Bool.(mask))
    keys = vec(all_indices)[vec(valid)]
    n = length(keys)
    n > 0 || throw(ArgumentError("No valid candidates after masking"))

    coords = _get_coordinates_from_cartesian_indices(keys)
    incl = _extract_and_process_inclusion(inclusion, keys, n)
    feat = Matrix(mat[valid]')

    return CandidatePool(
        n, 
        coords, 
        feat,
        incl, 
        keys,
        Dict{Symbol,Any}(:source => :matrix, :size => size(mat))
    )
end


# ========================================================================
# Tests
# ========================================================================

@testitem "We can construct CandidatePool from a matrix" begin
    mat = rand(10, 15)
    cp = CandidatePool(mat)
    @test cp.n == 150
    @test size(cp.coordinates) == (2, 150)
    @test cp.features == Matrix(vec(mat)')
    @test sum(cp.inclusion) ≈ 1.0
end

@testitem "We can construct CandidatePool from a matrix with mask" begin
    mask = falses(10, 10)
    mask[1:5, :] .= true
    cp = CandidatePool(rand(10, 10); mask)
    @test cp.n == 50
    @test all(k -> k[1] <= 5, cp.keys)
end

@testitem "We can construct CandidatePool from a matrix with custom inclusion" begin
    inclusion = rand(10, 10)
    cp = CandidatePool(rand(10, 10); inclusion)
    @test cp.n == 100
    @test cp.inclusion == vec(inclusion) ./ sum(inclusion)
end

@testitem "We can construct CandidatePool from a matrix with custom inclusion and a mask" begin
    inclusion = rand(10, 10)
    mask = falses(10, 10)
    mask[1:5, :] .= true
    cp = CandidatePool(rand(10, 10); inclusion, mask)
    @test cp.n == 50
    @test cp.inclusion == inclusion[mask] ./ sum(inclusion[mask])
end

@testitem "We can construct a CandidatePool from a BON" begin
    cp1 = CandidatePool(rand(10, 10))
    sr = SamplingResult(cp1, [1, 5, 10], SimpleRandom(n=3), Dict{Symbol,Any}())
    cp2 = candidatepool(sr)
    @test cp2.n == 3
    @test cp2.keys == sr.sites
end
