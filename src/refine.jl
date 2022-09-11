"""
    refine!(cooords::Vector{CartesianIndex}, pool::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})

Refines a set of candidate sampling locations in the preallocated vector `coords`
from a vector  of coordinates `pool` using `sampler`, where `sampler` is a [`BONRefiner`](@ref).
"""
function refine!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::ST,
    uncertainty::Matrix{T}
) where {ST <: BONRefiner, T <: AbstractFloat, N}
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
    return _generate!(coords, copy(pool), sampler, uncertainty)
end

"""
    refine!(cooords::Vector{CartesianIndex}, pool::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})

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
    return (p, u) -> refine!(coords, copy(p), sampler, u)
end

"""
    refine(pool::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})

Refines a set of candidate sampling locations and returns a vector `coords` of length numpoints
from a vector  of coordinates `pool` using `sampler`, where `sampler` is a [`BONRefiner`](@ref).
"""
function refine(
    pool::Vector{CartesianIndex},
    sampler::ST,
    uncertainty::Matrix{T},
) where {ST <: BONRefiner, T <: AbstractFloat}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return refine!(coords, copy(pool), sampler, uncertainty)
end

"""
    refine(sampler::BONRefiner)

Returns a curried function of `refine` with *two* methods: both are using the
output of `seed`, one in its packed form, the other in its splatted form.
"""
function refine(sampler::ST) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    _inner(p, u) = refine!(coords, copy(p), sampler, u)
    _inner(p) = refine!(coords, first(p), sampler, last(p))
    return _inner
end

"""
    refine(pack, sampler::BONRefiner)
    
Calls `refine` on the appropriatedly splatted version of `pack`.
"""
function refine(
    pack::Tuple{Vector{CartesianIndex}, Matrix{Float64}},
    sampler::ST,
) where {ST <: BONRefiner}
    return refine(first(pack), sampler, last(pack))
end