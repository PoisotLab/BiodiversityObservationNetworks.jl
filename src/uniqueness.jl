"""
    Uniqueness

A `BONRefiner`
"""
Base.@kwdef struct Uniqueness{I <: Integer, T <: AbstractFloat} <: BONRefiner
    numpoints::I = 30
    layers::Array{T, 3} = rand(50, 50, 3)
    function Uniqueness(numpoints, layers::Array{T, N}) where {T, N}
        uniq = new{typeof(numpoints), T}(numpoints, layers)
        check_arguments(uniq)
        return uniq
    end
end

function check_arguments(uniq::Uniqueness)
    check(TooFewSites, uniq)
    
    length(size(uniq.layers)) == 3 || 
        throw(
            ArgumentError(
                "You cannot have a Uniqueness sampler without layers passed as a 3-dimenisional array, where the first two axes are latitude and longitude.",
            ),
        )

    max_num_points = prod(size(uniq.layers)[1:2])
    max_num_points >= uniq.numpoints || throw(TooManySites("Number of requested sites $(uniq.numpoints) is more than the number of candidates $max_num_points."))
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::Uniqueness,
    uncertainty,
) where {T <: AbstractFloat}
    layers = sampler.layers
    ndims(layers) <= 2 &&
        throw(ArgumentError("Uniqueness needs more than one layer to work."))
    size(uncertainty) != (size(layers, 1), size(layers, 2)) &&
        throw(DimensionMismatch("Layers are not the same dimension as uncertainty"))

    covscore = zeros(length(pool))
    for (i, p1) in enumerate(pool)
        v1 = layers[p1[1], p1[2], :]
        for (j, p2) in enumerate(pool)
            v2 = layers[p2[1], p2[2], :]
            if p1 != p2
                covscore[i] += abs(cov(v1, v2))
            end
        end
    end

    np = sampler.numpoints
    sortedvals = sortperm(vec(covscore))

    coords[:] .= pool[sortedvals[1:np]]
    return (coords, uncertainty)
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
    @test_throws TooFewSites Uniqueness(numpoints=-1)
    @test_throws TooFewSites Uniqueness(numpoints=0)
    @test_throws TooFewSites Uniqueness(numpoints=1)
end

@testitem "Uniqueness throws error if more points are requested than are possible" begin
    @test_throws TooManySites Uniqueness(numpoints=26, layers=rand(5,5,2))
end

@testitem "Uniqueness has correct subtypes" begin
    @test typeof(Uniqueness(2, zeros(5, 5, 5))) <: BONRefiner
    @test typeof(Uniqueness(2, zeros(5, 5, 5))) <: BONSampler
end


