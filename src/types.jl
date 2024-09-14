"""
    abstract type BONSampler end

A `BONSampler` is any algorithm for proposing a set of sampling locations.
"""
abstract type BONSampler end

numsites(s::BONSampler) = s.numsites
pool(s::BONSampler) = s.pool

"""
    mutable struct Sites{T}

"""
mutable struct Sites{T}
    coordinates::Vector{T}
end
_allocate_sites(n) = Sites(Array{CartesianIndex}(undef, n))
coordinates(s::Sites) = s.coordinates
Base.getindex(s::Sites, i::Integer) = getindex(coordinates(s), i)
Base.setindex!(s::Sites, v, i::Integer) = setindex!(coordinates(s), v,i)
Base.length(s::Sites) = length(coordinates(s))
Base.eachindex(s::Sites) = eachindex(s.coordinates)

abstract type LayerType end 
abstract type DataLayer <: LayerType end 
abstract type InclusionProbability <: LayerType end 

struct Layer{T<:LayerType,L}
    layer::L
end 
pool(l::Layer) = Sites(vec(findall(l.layer.indices)))
Base.size(l::Layer) = size(l.layer)

struct Stack{T<:LayerType,N,L}
    layers::Dict{N,Layer{T,L}}
end
