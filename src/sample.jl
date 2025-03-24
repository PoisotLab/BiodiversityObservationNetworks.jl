"""
    BONSampler

An abstract type that is the supertype for all methods for sampling
BiodiversityObservationNetworks. 
"""
abstract type BONSampler end 

struct MultistageSampler <: BONSampler  
    samplers::Vector{<:BONSampler} 
end

const __BON_DOMAINS = Union{Vector{<:SDMLayer}, SDMLayer, Polygon, Vector{<:Polygon}, BiodiversityObservationNetwork}


Base.getindex(samplers::MultistageSampler, i::Integer) = samplers.samplers[i]
Base.firstindex(samplers::MultistageSampler) = firstindex(samplers.samplers)
Base.eachindex(samplers::MultistageSampler) = eachindex(samplers.samplers)
Base.iterate(samplers ::MultistageSampler) = iterate(samplers, firstindex(samplers.samplers))
Base.iterate(samplers ::MultistageSampler, i) = Base.iterate(samplers.samplers, i)


# so its possible that there will be a polygon and a raster/rasterstack 
# anything that works on a rasterstack can also work on a BON with set of
# covariates.


function sample(samplers::MultistageSampler, domain::__BON_DOMAINS)
    bon = sample(first(samplers), domain)
    for i in eachindex(samplers)[2:end]
        bon = sample(samplers[i], domain, bon)
    end
    return bon
end 

function _what_did_you_pass(geom)
    is_polygonizable(geom) && return Polygon
    #is_rasterizable(geom) && return Raster
    is_bonifyable(geom) && return BiodiversityObservationNetwork
    return nothing
end


"""
    sample

Sample from `geom`. This is highest level dispatch which assumes nothing about
the geometry the user is trying to apply [`BONSampler`](@ref) to.
"""
function sample(sampler::S, geom::T) where {S<:BONSampler,T}
    GEOM_TYPE = _what_did_you_pass(geom)
    isnothing(GEOM_TYPE) && throw(ArgumentError("$T cannot be coerced to a valid Geometry. Valid geometries for $S are $(_valid_geometries(sampler))"))
    sample(sampler, Base.convert(GEOM_TYPE, geom))
end

"""
    sample

Attempt to use `BONSampler` to sample from a valid `geom`
"""
sample(sampler::BONSampler, geom::__BON_DOMAINS) = _sample(sampler, geom)

"""
    sample

Attempt to use `BONSampler` to sample from a valid `geom`
"""
function sample(sampler::BONSampler, geom::__BON_DOMAINS, bon::BiodiversityObservationNetwork)
    _sample(sampler, geom, bon)
end
