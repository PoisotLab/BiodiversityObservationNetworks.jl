function seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T}) where {ST<:BONSampler,T<:AbstractFloat}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    return _generate!(coords, sampler, uncertainty)
end

function seed!(coords::Vector{CartesianIndex}, sampler::ST) where {ST<:BONSampler}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    return (u) -> _generate!(coords, sampler, u)
end

function seed(sampler::ST, uncertainty::Matrix{T}) where {ST <: BONSampler, T<:AbstractFloat}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    seed!(coords, sampler, uncertainty)
    return coords
end

function seed(sampler::ST) where {ST <: BONSampler}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return (u) -> _generate!(coords, sampler, u)
end
