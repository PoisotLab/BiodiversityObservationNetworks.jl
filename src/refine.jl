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
    if length(coords) != sampler.numsites
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numsites` fields of the sampler",
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
    if length(coords) != sampler.numsites
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numsites` fields of the sampler",
            ),
        )
    end
    return (p) -> refine!(coords, copy(p), sampler)
end

"""
    refine(pool::Vector{CartesianIndex}, sampler::ST)

Refines a set of candidate sampling locations and returns a vector `coords` of length numsites
from a vector  of coordinates `pool` using `sampler`, where `sampler` is a [`BONRefiner`](@ref).
"""
function refine(
    pool::Vector{CartesianIndex},
    sampler::ST,
) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numsites)
    return refine!(coords, copy(pool), sampler)
end

"""
    refine(sampler::BONRefiner)

Returns a curried function of `refine`
"""
function refine(sampler::ST) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numsites)
    _inner(p) = refine!(coords, first(p), sampler)
    return _inner
end

"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})

Puts a set of candidate sampling locations in the preallocated vector `coords`
from a raster `uncertainty` using `sampler`, where `sampler` is a [`BONSeeder`](@ref).

  - Seeder's work on rasters, refiners work on set of coordinates.
"""
function seed!(
    coords::Vector{CartesianIndex},
    sampler::ST,
) where {ST <: BONSeeder}
    length(coords) != sampler.numsites &&
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numsites` fiel s of the sampler",
            ),
        )
    return _generate!(coords, sampler)
end

"""
    seed(sampler::ST)

Produces a set of candidate sampling locations in a vector `coords` of length numsites
from a raster using `sampler`, where `sampler` is a [`BONSeeder`](@ref).
"""
function seed(sampler::ST) where {ST <: BONSeeder}
    coords = Vector{CartesianIndex}(undef, sampler.numsites)
    return seed!(coords, sampler)
end
