"""
    voronoi(bon, domain)

Construct Voronoi polygons for the nodes in a
`BiodiversityObservationNetwork` within the given domain. The domain is coerced
to a `PolygonDomain` when needed. Output polygons are clipped to the domain
extent.
"""
voronoi(bon, domain) = voronoi(bon, to_polygon(domain))
voronoi(bon, domain::RasterDomain)  = voronoi(bon, to_polygon(domain))

voronoi(bon, domain::SimpleSDMPolygons.AbstractGeometry) = voronoi(bon, to_polygon(domain))

function voronoi(bon, domain::PolygonDomain)
    (xmin, xmax), (ymin, ymax) = extent(domain)
    bbox = (xmin, xmax, ymin, ymax) 

    coord = getcoordinates(bon)

    tri = DT.triangulate(coord)
    vor = DT.voronoi(tri)

    polys = [intersect(domain.data, LeanBONsAPI.Polygon(DT.get_polygon_coordinates(vor, i, bbox))) for i in DT.each_polygon_index(vor)]
    return polys
end