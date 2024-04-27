"""
    Uniqueness

A `BONRefiner`
"""
Base.@kwdef mutable struct Uniqueness{I <: Integer, T <: Number} <: BONRefiner
    numpoints::I = 30
    layers::Array{T, 3}
    function Uniqueness(numpoints, layers::Array{T, N}) where {T, N}
        if numpoints < one(numpoints)
            throw(
                ArgumentError(
                    "You cannot have a Uniqueness sampler with less than one point",
                ),
            )
        end
        if N != 3
            throw(
                ArgumentError(
                    "You cannot have a Uniqueness sampler without layers passed as a cube.",
                ),
            )
        end
        return new{typeof(numpoints), T}(numpoints, layers)
    end
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
