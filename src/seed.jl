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
    length(coords) != sampler.numpoints &&
        throw(
            DimensionMismatch(
                "The length of the coordinate vector must match the `numpoints` fiel s of the sampler",
            ),
        )
    return _generate!(coords, sampler)
end

"""
    seed(sampler::ST)

Produces a set of candidate sampling locations in a vector `coords` of length numpoints
from a raster using `sampler`, where `sampler` is a [`BONSeeder`](@ref).
"""
function seed(sampler::ST) where {ST <: BONSeeder}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return seed!(coords, sampler)
end
