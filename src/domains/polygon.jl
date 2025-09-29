struct PolygonDomain{T}
    data::T
end

"""
    extent(p::PolygonDomain{T<:SimpleSDMPolygons.AbstractGeometry})
"""
function extent(p::PolygonDomain{<:SimpleSDMPolygons.AbstractGeometry})
    bbox = SimpleSDMPolygons.boundingbox(p.data)
    return Extent(X=(bbox.left, bbox.right), Y=(bbox.bottom, bbox.top))
end
extent(::PolygonDomain{T}) where {T} = error("extent is not implemented for Polygons of type $T")
extent(polys::Vector{<:PolygonDomain}) = union(extent.(polys)...)


contains(::PolygonDomain{T}, ::V) where {T,V} = error("contains is not implemented for Polygons of type $T and points of type $V")
contains(doms::Vector{<:PolygonDomain}, pt) = any(d -> contains(d, pt), doms)

Base.in(pt, domain::PolygonDomain) = contains(domain, pt)
Base.in(pt, doms::Vector{<:PolygonDomain}) = contains(doms, pt)

convert_domain(geom::SimpleSDMPolygons.AbstractGeometry) = PolygonDomain(geom)
convert_domain(domain::PolygonDomain) = domain


function grid_polygon(polygon, grid_size)
    x,y = extent(polygon)
    layer = SDMLayer(ones(Bool, grid_size), x=x, y=y)
    RasterDomain(mask(layer, polygon.data))
end

function grid_polygon(raster::RasterDomain, polygon)
    x,y = extent(raster)
    layer = SDMLayer(ones(Bool, size(raster)), x=x, y=y)
    RasterDomain(mask(layer, polygon.data))
end