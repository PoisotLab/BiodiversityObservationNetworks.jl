

function voronoi(bon::BiodiversityObservationNetwork, geom)
    tri = DT.triangulate(GI.coordinates(bon))
    vor = DT.voronoi(tri)

    (xmin, xmax), (ymin, ymax) = GI.extent(geom) 
    bbox = (xmin, xmax, ymin, ymax) 
    

    vor_polys = convert.(Polygon, [AGDAL.intersection(geom.geometry.geom, AGDAL.createpolygon(DT.get_polygon_coordinates(vor, i, bbox))) for i in DT.each_polygon_index(vor)])

    # This is O(n²) and I hate it
    # But necessary to make sure the Voronoi polygons are in the same order as the BON nodes
    order = [findfirst(p->GO.contains(v,p), [n.coordinate for n in bon]) for v in vor_polys]

    return vor_polys[order]
end 
