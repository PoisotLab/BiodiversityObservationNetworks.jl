"""
    gadm

A wrapper around the SpeciesDistributionToolkit `gadm` method. Ensures the
resulting object is using the ArchGDAL representation, which makes intersecting
Polygon's easier later on. 
"""
function gadm(arg)
    json_str = SDT.GeoJSON.write(SDT.gadm(arg)[1].geometry)
    poly = GI.MultiPolygon(AGDAL.fromJSON(json_str))
    Polygon(GI.trait(poly), poly)
end

function gadm(arg...)
    json_str = SDT.GeoJSON.write.(SDT.gadm(arg...))
    poly = GI.MultiPolygon.(AGDAL.fromJSON.(json_str))
    Polygon(GI.trait(first(poly)), poly)
end