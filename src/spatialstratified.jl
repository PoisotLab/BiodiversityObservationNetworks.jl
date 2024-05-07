"""
    SpatiallyStratified
"""
@kwdef struct SpatiallyStratified{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numsites::I = 50
    strata::Matrix{I} = _default_strata((50, 50))
    inclusion_probability_by_stratum::Vector{F} = ones(3) ./ 3
    function SpatiallyStratified(numsites, strata, inclusion_probability_by_stratum)
        ss = new{typeof(numsites), typeof(inclusion_probability_by_stratum[begin])}(
            numsites,
            strata,
            inclusion_probability_by_stratum,
        )
        check_arguments(ss)
        return ss
    end
end

maxsites(ss::SpatiallyStratified) = prod(size(ss.strata))

function check_arguments(ss::SpatiallyStratified)
    check(TooFewSites, ss)
    check(TooManySites, ss)

    length(unique(ss.strata)) == length(ss.inclusion_probability_by_stratum) || throw(
        ArgumentError(
            "Inclusion probability vector does not have the same number of strata as there are unique values in the strata matrix",
        ),
    )

    return sum(ss.inclusion_probability_by_stratum) ≈ 1.0 ||
           throw(ArgumentError("Inclusion probabilities for each strata do not sum to 1."))
end

function _default_strata(sz)
    mat = zeros(typeof(sz[1]), sz...)

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
)
    strata = sampler.strata
    idx_per_strata = [
        findall(i -> strata[i] == x, CartesianIndices(strata)) for
        x in unique(sampler.strata)
    ]
    πᵢ = sampler.inclusion_probability_by_stratum

    strata_per_sample = rand(Categorical(πᵢ), sampler.numsites)
    for (i, s) in enumerate(strata_per_sample)
        coords[i] = rand(idx_per_strata[s])
    end

    return coords
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
    coords = seed(ss)
    @test typeof(coords) <: Vector{CartesianIndex}
end

@testitem "SpatiallyStratified throws error when number of sites is below 2" begin
    @test_throws TooFewSites SpatiallyStratified(numsites = -1)
    @test_throws TooFewSites SpatiallyStratified(numsites = 0)
    @test_throws TooFewSites SpatiallyStratified(numsites = 1)
end

@testitem "SpatiallyStratified can use custom number of points as keyword argument" begin
    NUM_POINTS = 42
    ss = SpatiallyStratified(; numsites = NUM_POINTS)
    @test ss.numsites == NUM_POINTS
end

@testitem "SpatiallyStratified can use custom strata as keyword argument" begin
    dims = (42, 30)
    strata = rand(1:10, dims...)
    inclusion_probability = [0.1 for i in 1:10]
    ss = SpatiallyStratified(;
        strata = strata,
        inclusion_probability_by_stratum = inclusion_probability,
    )
    coords = seed(ss)
    @test typeof(coords) <: Vector{CartesianIndex}
end

@testitem "SpatiallyStratified throws error if the number of inclusion probabilities are different than the number of unique strata" begin
    dims = (42, 42)
    inclusion_probability = [0.5, 0.5]
    strata = rand(1:5, dims...)

    @test_throws ArgumentError SpatiallyStratified(;
        strata = strata,
        inclusion_probability_by_stratum = inclusion_probability,
    )
end
