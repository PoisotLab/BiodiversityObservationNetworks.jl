#=

Mask conversion and application utilities.

`convert_mask` coerces mask inputs (matrices, `SDMLayer`s, polygons, or
`RasterDomain`s) to a mask aligned with a given raster `domain`, performing
consistency checks on size, extent, and CRS where applicable. `mask!` applies
the mask in-place to a raster domain or each raster in a stack by updating
`domain.pool`/indices.

=#

# Nothing -> Nothing
convert_mask(::Any, ::Nothing; kwargs...) = nothing
convert_mask(::RasterStack, ::Nothing; kwargs...) = nothing

convert_mask(rs::RasterStack, mask; kwargs...) = convert_mask(first(rs.rasters), mask; kwargs...) 


"""
    convert_mask(domain::RasterDomain, mask::AbstractMatrix; kwargs...)

Convert a boolean matrix mask to a `RasterDomain`-aligned mask.
"""
function convert_mask(domain::RasterDomain, mask::AbstractMatrix; kwargs...)
    # TODO: assert matrix is bool
    convert_mask(domain, RasterDomain(ones(size(mask)), Matrix(mask)))
end 

"""
    convert_mask(domain::RasterDomain, mask::SDMLayer; kwargs...)

Convert an `SDMLayer` to a `RasterDomain` mask with aligned indices.
"""
function convert_mask(domain::RasterDomain, mask::SDMLayer; kwargs...)
    convert_mask(domain, RasterDomain(mask, mask.indices); kwargs...)
end 

"""
    convert_mask(domain::RasterDomain{<:SDMLayer}, mask::RasterDomain{<:SDMLayer}; kwargs...)

Validate size, extent, and CRS when both domain and mask are SDMLayer-backed.
"""
function convert_mask(domain::RasterDomain{<:SDMLayer}, mask::RasterDomain{<:SDMLayer}; kwargs...)
    size(domain) == size(mask) || throw(ArgumentError("Domain and Mask must be same size."))
    extent(domain) == extent(mask) || throw(ArgumentError("Domain and Mask must have same extent."))
    crs(domain) == crs(mask) || throw(ArgumentError("Domain and Mask must have CRS."))

    return mask
end 

"""
    convert_mask(domain::RasterDomain{T}, mask::RasterDomain{<:Matrix}; kwargs...) where T

Validate size when the mask is matrix-backed.
"""
function convert_mask(domain::RasterDomain{T}, mask::RasterDomain{<:Matrix}; kwargs...) where T<:Union{<:SDMLayer, <:Matrix}
    size(domain) == size(mask) || throw(ArgumentError("Domain and Mask must be same size."))
    return mask
end 

"""
    convert_mask(domain::RasterDomain, mask::SimpleSDMPolygons.AbstractGeometry)

Coerce a polygon geometry to a `RasterDomain` mask by rasterizing the polygon.
"""
function convert_mask(domain::RasterDomain, mask::SimpleSDMPolygons.AbstractGeometry)
    return convert_mask(domain, PolygonDomain(mask))
end

"""
    convert_mask(domain::RasterDomain, mask::PolygonDomain; kwargs...)

Rasterize a polygon mask into the grid of `domain` and return a
`RasterDomain`-aligned mask.
"""
function convert_mask(domain::RasterDomain, mask::PolygonDomain; kwargs...)
    # Create full SDMLayer
    grid = zeros(Float64, size(domain))
    layer = RasterDomain(SDMLayer(grid; x=domain.data.x, y=domain.data.y, crs=domain.data.crs))
    
    # Mask it by polygon
    mask!(layer, mask)
    return layer
end


"""
    mask!(domain, mask)

Apply `mask` to `domain` in-place. Overloads handle `RasterDomain`s backed by
`SDMLayer` or `Matrix`, and propagate to all rasters in a `RasterStack`.
Polygon masks are supported for SDMLayer-backed rasters.
"""
mask!(::Any, ::Nothing) = nothing
mask!(::RasterStack, ::Nothing) = nothing

function mask!(rs::RasterStack, args...)
    for r in rs.rasters
        mask!(r, args...)
    end
end

function mask!(::RasterDomain{<:Matrix}, ::Union{SimpleSDMPolygons.AbstractGeometry,PolygonDomain})
    throw(ArgumentError("Matrix based RasterDomains cannot be masked with polygons"))
end

function mask!(domain::RasterDomain{<:SDMLayer}, _mask::SimpleSDMPolygons.AbstractGeometry)
    SimpleSDMLayers.mask!(domain.data, _mask)
    domain.pool .= domain.data.indices 
end 

function mask!(domain::RasterDomain{<:SDMLayer}, _mask::PolygonDomain)
    SimpleSDMLayers.mask!(domain.data, _mask.data)
    domain.pool .= domain.data.indices 
end 


function mask!(
    domain::RasterDomain{<:Matrix}, 
    _mask::RasterDomain{<:Matrix}
)
    domain.pool .= _mask.pool 
end 


function mask!(
    domain::RasterDomain{<:SDMLayer}, 
    _mask::RasterDomain{<:Matrix}
)
    domain.data.indices .= domain.data.indices .&& _mask.pool
    domain.pool .= domain.data.indices
end 


function mask!(
    domain::RasterDomain{<:SDMLayer}, 
    _mask::RasterDomain{<:SDMLayer}
)
    domain.data.indices .= _mask.data.indices
    SimpleSDMLayers.mask!(domain.data, _mask.data)
    domain.pool .= domain.data.indices
end 


# Dealing with BONs
function mask!(
    domain::BiodiversityObservationNetwork, 
    _mask::RasterDomain{<:SDMLayer}
)
    
end 

function mask!(
    domain::BiodiversityObservationNetwork, 
    _mask::SimpleSDMPolygons.AbstractGeometry
)
    
end 
