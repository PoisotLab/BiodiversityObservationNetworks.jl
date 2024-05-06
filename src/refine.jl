"""
    refine!(cooords::Vector{CartesianIndex}, pool::Vector{CartesianIndex}, sampler::ST)

Refines a set of candidate sampling locations in the preallocated vector `coords`
from a vector of coordinates `pool` using `sampler`, where `sampler` is a [`BONRefiner`](@ref).
"""
function refine!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::ST,
) where {ST <: BONRefiner}
    if length(coords) != sampler.numpoints
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numpoints` fields of the sampler",
            ),
        )
    end
    if length(coords) > length(pool)
        throw(
            DimensionMismatch(
                "The number of refined points must be at least the number of seeded points",
            ),
        )
    end
    return _generate!(coords, copy(pool), sampler)
end

"""
    refine!(cooords::Vector{CartesianIndex}, pool::Vector{CartesianIndex}, sampler::ST)

The curried version of `refine!`, which returns a function that acts on the input
coordinate pool passed to the curried function (`p` below).
"""
function refine!(coords::Vector{CartesianIndex}, sampler::ST) where {ST <: BONRefiner}
    if length(coords) != sampler.numpoints
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numpoints` fields of the sampler",
            ),
        )
    end
    return (p) -> refine!(coords, copy(p), sampler)
end

"""
    refine(pool::Vector{CartesianIndex}, sampler::ST)

Refines a set of candidate sampling locations and returns a vector `coords` of length numpoints
from a vector  of coordinates `pool` using `sampler`, where `sampler` is a [`BONRefiner`](@ref).
"""
function refine(
    pool::Vector{CartesianIndex},
    sampler::ST,
) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return refine!(coords, copy(pool), sampler)
end

"""
    refine(sampler::BONRefiner)

Returns a curried function of `refine`
"""
function refine(sampler::ST) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    _inner(p) = refine!(coords, first(p), sampler)
    return _inner
end
