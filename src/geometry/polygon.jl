struct Polygon{T, G}
    geometry::G
    function Polygon(::T, geom::G) where {T <: GI.AbstractTrait, G}
        return new{T, G}(geom)
    end
end
Base.show(io, ::Polygon) = print(io, "Polygon")
Base.show(io, vec::Vector{<:Polygon}) = print(io, "Vector of $(length(vec)) Polygons")


GI.isgeometry(::Polygon)::Bool = true
GI.geomtrait(::Polygon)::DataType = GI.MultiPolygonTrait
GI.ngeom(::GI.MultiPolygonTrait, geom::Polygon)::Integer = ngeom(geom.geometry)
GI.getgeom(::GI.MultiPolygonTrait, geom::Polygon, i) = getgeom(GI.MultiPolygonTrait, geom.geometry, i)
GI.extent(geom::Polygon) = GI.extent(geom.geometry)

GO.contains(geom::Polygon, coord) = GO.contains(geom.geometry, coord)
GO.area(geom::Polygon) = GO.area(geom.geometry)

Base.convert(::Type{Polygon}, foo) = _convert_to_bons_polygon(foo)

const __POLYGONIZABLE_TYPES = Union{<:GJSON.FeatureCollection,<:GJSON.MultiPolygon,Vector{<:GJSON.MultiPolygon}}

is_polygonizable(::T) where T = T <: __POLYGONIZABLE_TYPES


# =================================================================
# GeoJSON handlers
# 
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
