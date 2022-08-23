"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T})


Puts a set of candidate sampling locations in the preallocated vector `coords` 
from a raster `uncertainty` using `sampler`, where `sampler` is a [`BONSeeder`](@ref).


- Seeder's work on rasters, refiners work on set of coordinates. 
"""
function seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T}) where {ST<:BONSampler,T<:AbstractFloat}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    return _generate!(coords, sampler, uncertainty)
end

"""
    seed!(coords::Vector{CartesianIndex}, sampler::ST)

The curried version of `seed!`, which returns a function that acts on the input
uncertainty layer passed to the curried function (`u` below).
"""
function seed!(coords::Vector{CartesianIndex}, sampler::ST) where {ST<:BONSampler}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    return (u) -> _generate!(coords, sampler, u)
end

"""
    seed(sampler::ST, uncertainty::Matrix{T})


"""
function seed(sampler::ST, uncertainty::Matrix{T}) where {ST <: BONSampler, T<:AbstractFloat}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    seed!(coords, sampler, uncertainty)
    return coords
end

function seed(sampler::ST) where {ST <: BONSampler}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return (u) -> _generate!(coords, sampler, u)
end
