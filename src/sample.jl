abstract type BONSampler end 

function _what_did_you_pass(geom)
    is_polygonizable(geom) && return Polygon
    is_rasterizable(geom) && return Raster
    is_bonifyable(geom) && return BiodiversityObservationNetwork
end

const __BON_DOMAINS = Union{Raster, RasterStack, Polygon, Vector{<:Polygon}, BiodiversityObservationNetwork}

function sample(sampler::BONSampler, geom::Any)
    T = _what_did_you_pass(geom)
    sample(sampler, Base.convert(T, geom))
end

function sample(sampler::BONSampler, geom::__BON_DOMAINS)
    _sample(sampler, geom)
end
