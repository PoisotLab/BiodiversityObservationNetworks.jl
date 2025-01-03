"""
    BONSampler

An abstract type that is the supertype for all methods for sampling
BiodiversityObservationNetworks. 
"""
abstract type BONSampler end 

function _what_did_you_pass(geom)
    is_polygonizable(geom) && return Polygon
    is_rasterizable(geom) && return Raster
    is_bonifyable(geom) && return BiodiversityObservationNetwork
end

const __BON_DOMAINS = Union{Raster, RasterStack, Polygon, Vector{<:Polygon}, BiodiversityObservationNetwork}

"""
    sample

Sample from `geom`. This is highest level dispatch which assumes nothing about
the geometry the user is trying to apply [`BONSampler`](@ref) to.
"""
function sample(sampler::BONSampler, geom::T) where T
    GEOM_TYPE = _what_did_you_pass(geom)
    isnothing(GEOM_TYPE) && throw(ArgumentError("$T cannot be coerced to a valid Geometry"))
    sample(sampler, Base.convert(GEOM_TYPE, geom))
end

"""
    sample

Attempt to use `BONSampler` to sample from a valid `geom`
"""
function sample(sampler::BONSampler, geom::__BON_DOMAINS)
    _sample(sampler, geom)
end
