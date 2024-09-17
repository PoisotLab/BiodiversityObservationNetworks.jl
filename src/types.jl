"""
    abstract type BONSampler end

A `BONSampler` is any algorithm for proposing a set of sampling locations.
"""
abstract type BONSampler end

numsites(s::BONSampler) = s.numsites

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
struct CategoricalData <: DataLayer end 
struct ContinousData <: DataLayer end 

abstract type InclusionProbability <: LayerType end 
abstract type ResultLayer <: LayerType end 


# Distribution over categorical values at each pixel
# Scalar at each pixel 
struct CategoricalResult <: LayerType end 
struct ContinuousResult <: LayerType end 

# Layer
struct Layer{T<:LayerType,L}
    layer::L
end 
Base.size(l::Layer) = size(l.layer)
function Layer(l::S; categorical = false) where S<:SimpleSDMLayers.SDMLayer
    T = categorical ? CategoricalData : ContinousData
    Layer{T,S}(l)
end 
Layer(m::Matrix{I}) where I<:Integer = Layer{CategoricalData,typeof(m)}(m)

_indices(m::M) where M<:Matrix = findall(i->!isnothing(m[i]) && !isnan(m[i]) && !ismissing(m[i]), CartesianIndices(m))
_indices(l::SDMLayer) = findall(l.indices)

# Pool method 
pool(l::Layer) = Sites(vec(_indices(l.layer)))
pool(dims::Tuple{I,I}) where I<:Integer = Sites(vec(collect(CartesianIndices(dims))))


# Layer set w/ names
struct Stack{T<:LayerType,N,L}
    layers::Dict{N,Layer{T,L}}
end
