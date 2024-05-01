"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})

Puts a set of candidate sampling locations in the preallocated vector `coords`
from a raster `uncertainty` using `sampler`, where `sampler` is a [`BONSeeder`](@ref).

  - Seeder's work on rasters, refiners work on set of coordinates.
"""
function seed!(
    coords::Vector{CartesianIndex},
    sampler::ST,
    uncertainty::Matrix{T},
) where {ST <: BONSeeder, T <: AbstractFloat}
    length(coords) != sampler.numpoints &&
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numpoints` fiel s of the sampler",
            ),
        )

    max_sites = prod(size(uncertainty))
    max_sites < sampler.numpoints && throw(
        TooManySites(
            "Cannot select $(sampler.numpoints) sampling sites from $max_sites possible locations.",
        ),
    )
    return _generate!(coords, sampler, uncertainty)
end

"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST)

The curried version of `seed!`, which returns a function that acts on the input
uncertainty layer passed to the curried function (`u` below).
"""
function seed!(coords::Vector{CartesianIndex}, sampler::ST) where {ST <: BONSeeder}
    if length(coords) != sampler.numpoints
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numpoints` fields of the sampler",
            ),
        )
    end
    return (u) -> seed!(coords, sampler, u)
end

"""
    seed(sampler::ST, uncertainty::Matrix{T})

Produces a set of candidate sampling locations in a vector `coords` of length numpoints
from a raster `uncertainty` using `sampler`, where `sampler` is a [`BONSeeder`](@ref).
"""
function seed(
    sampler::ST,
    uncertainty::Matrix{T},
) where {ST <: BONSeeder, T <: AbstractFloat}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return seed!(coords, sampler, uncertainty)
end

"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST)

The curried version of `seed!`, which returns a function that acts on the input
uncertainty layer passed to the curried function (`u` below).
"""
function seed(sampler::ST) where {ST <: BONSeeder}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return (u) -> seed!(coords, sampler, u)
end
