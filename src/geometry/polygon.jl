
abstract type Geometry end

"""
    Polygon

A Polygon extends the GeoInterface API for both PolygonTrait and
MultiPolygonTraits. 
"""
struct Polygon <: Geometry
    geometry::AG.IGeometry{AG.wkbMultiPolygon}
    function Polygon(poly::AG.IGeometry{AG.wkbMultiPolygon}) 
        new(poly)
    end 
end
const __POLYGONIZABLE_TYPES = Union{<:GJ.FeatureCollection,<:GJ.MultiPolygon,Vector{<:GJ.MultiPolygon}, <:AG.IGeometry{AG.wkbPolygon}, <:GI.Wrappers.Polygon}

# Constructors
Polygon(multipolys::Vector{<:AG.IGeometry{AG.wkbMultiPolygon}}) = Polygon.(multipolys)
Polygon(poly::AG.IGeometry{AG.wkbPolygon}) = Polygon(AG.forceto(poly, AG.wkbMultiPolygon))
Polygon(polys::Vector{<:AG.IGeometry{AG.wkbPolygon}}) = Polygon.(polys)
Polygon(feat::GJ.FeatureCollection) = Polygon(AG.fromJSON(GJ.write(feat.geometry[1])))
Polygon(multipoly::GJ.MultiPolygon) = Polygon(AG.fromJSON(GJ.write(multipoly)))
Polygon(res::Vector{<:GJ.MultiPolygon}) = Polygon.(AG.fromJSON.(GJ.write.(res)))

Base.convert(Polygon, t::__POLYGONIZABLE_TYPES) = Polygon(t)

is_polygonizable(::T) where T = T <: __POLYGONIZABLE_TYPES

# GeoInterface overloads
GI.isgeometry(::Polygon)::Bool = true
GI.geomtrait(::Polygon)::DataType = GI.MultiPolygonTrait
GI.ngeom(::Type{GI.MultiPolygonTrait}, geom::Polygon)::Integer = GI.ngeom(geom.geometry)
GI.getgeom(::Type{GI.MultiPolygonTrait}, geom::Polygon, i) = 
GI.getgeom(geom.geometry, i)
GI.crs(::Type{GI.MultiPolygonTrait}, geom::Polygon)= GI.crs(geom.geometry)
GI.extent(::Type{GI.MultiPolygonTrait}, geom::Polygon)::GI.Extents.Extent = GI.extent(geom.geometry)

# GeometryOps.jl overloads
GO.contains(geom::Polygon, pt) = GO.contains(geom.geometry, pt) 

# SDT Overloads
SDT.boundingbox(poly::Polygon) = begin
    (xm,xM), (ym, yM) = GI.extent(poly)
    return (left=xm, right=xM, bottom=ym, top=yM)
end 

SDT.mask!(layers, poly::Polygon) = begin
    geojson_poly = GJ.read(AG.toJSON(poly.geometry))
    SDT.mask!(layers, geojson_poly)
end

#=
struct Polygon{T, G}
    geometry::G
    function Polygon(::T, geom::G) where {T <: GI.AbstractTrait, G}
        return new{T, G}(geom)
    end
end
Base.show(io::IO, ::Polygon{T,V}) where {T,V} = print(io, "Polygon with inner geometry $V")
Base.show(io::IO, vec::Vector{<:Polygon}) = print(io, "Vector of $(length(vec)) Polygons")


GI.isgeometry(::Polygon)::Bool = true
GI.geomtrait(::Polygon)::DataType = GI.MultiPolygonTrait
GI.ngeom(::GI.MultiPolygonTrait, geom::Polygon)::Integer = ngeom(geom.geometry)
GI.getgeom(::GI.MultiPolygonTrait, geom::Polygon, i) = getgeom(GI.MultiPolygonTrait, geom.geometry, i)
GI.extent(geom::Polygon) = GI.extent(geom.geometry)

GO.contains(geom::Polygon, coord) = GO.contains(geom.geometry, coord)
GO.area(geom::Polygon) = GO.area(geom.geometry)


SDT.SimpleSDMLayers.mask!(layers, geom) = SDT.SimpleSDMLayers.mask!(layers, SDT.GeoJSON.read(AG.toJSON(geom.geometry.geom)))
SDT.boundingbox(geom::Polygon) = begin
    (xm,xM), (ym,yM)  = GI.extent(geom)
    return (left=xm, right=xM, bottom=ym, top=yM)
end 

Base.convert(::Type{Polygon}, foo) = _convert_to_bons_polygon(foo)

const __POLYGONIZABLE_TYPES = Union{<:GJSON.FeatureCollection,<:GJSON.MultiPolygon,Vector{<:GJSON.MultiPolygon}, <:AG.IGeometry{AG.wkbPolygon}, <:GI.Wrappers.Polygon}
is_polygonizable(::T) where T = T <: __POLYGONIZABLE_TYPES


# =================================================================
# GeoJSON handlers
# 
_convert_to_bons_polygon(geom::GI.Wrappers.Polygon) = Polygon(GI.trait(geom), geom)
function _convert_to_bons_polygon(geom::SDT.GeoJSON.MultiPolygon)
    Polygon(GI.trait(geom), geom)
end
function _convert_to_bons_polygon(vec::Vector{<:SDT.GeoJSON.MultiPolygon})
    if length(vec) > 1
        return [_convert_to_bons_polygon(vec[i]) for i in eachindex(vec)]
    else
        return _convert_to_bons_polygon(first(vec))
    end 
end 
function _convert_to_bons_polygon(fc::SDT.GeoJSON.FeatureCollection)
    try 
        return _convert_to_bons_polygon(fc.geometry)
    catch e 
        @error e
    end  
end 


# =================================================================
# ArchGDAL handlers
# 
function _convert_to_bons_polygon(geom::AG.IGeometry{AG.wkbPolygon})
    Polygon(GI.trait(geom), geom)
end

function _convert_to_bons_polygon(geom::AG.IGeometry{AG.wkbMultiPolygon})
    Polygon(GI.trait(geom), geom)
end
=#