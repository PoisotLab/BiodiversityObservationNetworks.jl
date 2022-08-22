@kwdef mutable struct BiodiversityObservationNetwork{IT<:Integer,M<:AbstractMatrix}
    numobservatories::IT = 50
    coordinates::M = missing
end

abstract type BONSeeder end
abstract type BONRefiner end

BONSampler = Union{BONSeeder,BONRefiner}

function _checknumpoints(sampler::T) where {T<:BONSampler}
    getfield(sampler, :numpoints) > 0 || throw(ArgumentError("numpoints boooooo"))
    return nothing
end


function seed!(coords::Vector{CartesianIndex}, sampler::ST, uncertainty::Matrix{T}) where {ST <: BONSampler, T<:AbstractFloat}
    return _generate!(coords, sampler, uncertainty)
end

function seed!(coords::Vector{CartesianIndex}, sampler::ST) where {ST <: BONSampler}
    return (u) -> _generate!(coords, sampler, u)
end

function seed(sampler::ST, uncertainty::Matrix{T}) where {ST <: BONSampler, T<:AbstractFloat}
    coords = CartesianIndex[] # TODO size
    seed!(coords, sampler, uncertainty)
    return coords
end

function seed(sampler::ST) where {ST <: BONSampler}
    coords = CartesianIndex[] # TODO size
    return (u) -> _generate!(coords, sampler, u)
end
