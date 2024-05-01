"""
    SpatiallyStratified
"""
@kwdef struct SpatiallyStratified{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numpoints::I = 50
    strata::Matrix{I} = _default_strata((50, 50))
    inclusion_probability_by_stratum::Vector{F} = ones(3) ./ 3
end

function _default_strata(sz)
    mat = zeros(Int64, sz...)

    x = sz[1] ÷ 2
    y = sz[2] ÷ 3

    mat[begin:x, :] .= 1
    mat[(x + 1):end, begin:y] .= 2
    mat[(x + 1):end, (y + 1):end] .= 3

    return mat
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::SpatiallyStratified,
    uncertainty::Matrix{T},
) where {T}
    strata = sampler.strata
    idx_per_strata = [
        findall(i -> strata[i] == x, CartesianIndices(strata)) for
        x in unique(sampler.strata)
    ]
    πᵢ = sampler.inclusion_probability_by_stratum

    strata_per_sample = rand(Categorical(πᵢ), sampler.numpoints)
    for (i, s) in enumerate(strata_per_sample)
        coords[i] = rand(idx_per_strata[s])
    end

    return coords, uncertainty
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "SpatiallyStratified default constructor works" begin
    @test typeof(SpatiallyStratified()) <: SpatiallyStratified
end

@testitem "SpatiallyStratified with default arguments can generate points" begin
    ss = SpatiallyStratified()
    uncert = rand(size(ss.strata)...)
    coords = seed(ss, uncert) |> first
    @test typeof(coords) <: Vector{CartesianIndex}
end
