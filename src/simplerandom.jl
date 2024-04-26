"""
    SimpleRandom

Implements Simple Random spatial sampling (each location has equal probability of selection)
"""
Base.@kwdef struct SimpleRandom{I <: Integer} <: BONSeeder
    numpoints::I = 50
    function SimpleRandom(numpoints)
        if numpoints < one(numpoints)
            throw(
                ArgumentError(
                    "You cannot have a SimpleRandom seeder with fewer than one point",
                ),
            )
        end
        return new{typeof(numpoints)}(numpoints)
    end
end

function _generate!(
    coords::Vector{CartesianIndex},
    sampler::SimpleRandom,
    uncertainty::Matrix{T},
) where {T <: AbstractFloat}
    pool = CartesianIndices(uncertainty)

    coords .= sample(pool, sampler.numpoints; replace = false)
    return (coords, uncertainty)
end
