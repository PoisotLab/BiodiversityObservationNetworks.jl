"""
    gadm

A wrapper around the SpeciesDistributionToolkit `gadm` method. Ensures the
resulting object is using the ArchGDAL representation, which makes intersecting
Polygon's easier later on. 
"""
function gadm(args...)
    Polygon(SDT.gadm(args...))
end
