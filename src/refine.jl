function refine!(coords::Vector{CartesianIndex}, pool::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T}) where {ST<:BONRefiner,T<:AbstractFloat}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    if length(coords) > length(pool)
        throw(DimensionMismatch("The number of refined points must be at least the number of seeded points"))
    end
    return _generate!(coords, copy(pool), sampler, uncertainty)
end

function refine!(coords::Vector{CartesianIndex}, sampler::ST) where {ST<:BONRefiner}
    if length(coords) != sampler.numpoints
        throw(DimensionMismatch("The length of the coordinate vector must match the `numpoints` fields of the sampler"))
    end
    if length(coords) > length(pool)
        throw(DimensionMismatch("The number of refined points must be at least the number of seeded points"))
    end
    return (p,u) -> _generate!(coords, p, sampler, u)
end

function refine(pool::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T}) where {ST <: BONRefiner, T<:AbstractFloat}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    refine!(coords, copy(pool), sampler, uncertainty)
    return coords
end

function refine(sampler::ST) where {ST <: BONRefiner}
    coords = Vector{CartesianIndex}(undef, sampler.numpoints)
    return (p, u) -> _generate!(coords, p, sampler, u)
end
