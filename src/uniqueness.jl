"""
    Uniqueness

A `BONSampler`
"""
Base.@kwdef struct Uniqueness{I <: Integer, T <: AbstractFloat} <: BONSampler
    numsites::I = 30
    layers::Array{T, 3} = rand(50, 50, 3)
    function Uniqueness(numsites, layers)
        uniq = new{typeof(numsites), typeof(layers[begin])}(numsites, layers)
        check_arguments(uniq)
        return uniq
    end
end

maxsites(uniq::Uniqueness) = prod(size(uniq.layers)[1:2])

function check_arguments(uniq::Uniqueness)
    check(TooFewSites, uniq)
    check(TooManySites, uniq)

    length(size(uniq.layers)) == 3 ||
        throw(
            ArgumentError(
                "You cannot have a Uniqueness sampler without layers passed as a 3-dimenisional array, where the first two axes are latitude and longitude.",
            ),
        )

    return nothing
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::Uniqueness,
)
    layers = sampler.layers
    ndims(layers) <= 2 &&
        throw(ArgumentError("Uniqueness needs more than one layer to work."))
    size(uncertainty) != (size(layers, 1), size(layers, 2)) &&
        throw(DimensionMismatch("Layers are not the same dimension as uncertainty"))

    total_covariance = zeros(length(pool))
    for (i, p1) in enumerate(pool)
        v1 = layers[p1[1], p1[2], :]
        for (j, p2) in enumerate(pool)
            v2 = layers[p2[1], p2[2], :]
            if p1 != p2
                total_covariance[i] += abs(cov(v1, v2))
            end
        end
    end

    np = sampler.numsites
    sortedvals = sortperm(vec(total_covariance))

    coords[:] .= pool[sortedvals[1:np]]
    return coords
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "Uniqueness default constructor works" begin
    @test typeof(Uniqueness()) <: Uniqueness
end

@testitem "Uniqueness requires more than one point" begin
    @test_throws TooFewSites Uniqueness(numsites = -1)
    @test_throws TooFewSites Uniqueness(numsites = 0)
    @test_throws TooFewSites Uniqueness(numsites = 1)
end

@testitem "Uniqueness throws error if more points are requested than are possible" begin
    @test_throws TooManySites Uniqueness(numsites = 26, layers = rand(5, 5, 2))
end

@testitem "Uniqueness works with positional constructor" begin
    @test typeof(Uniqueness(2, rand(5, 5, 5))) <: Uniqueness
end
