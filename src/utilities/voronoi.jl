function voronoi(bon::BiodiversityObservationNetwork, geom)
    tri = DT.triangulate(GI.coordinates(bon))
    vor = DT.voronoi(tri)

    (xmin, xmax), (ymin, ymax) = GI.extent(geom) 
    bbox = (xmin, xmax, ymin, ymax) 
    
    return convert.(Polygon, [AGDAL.createpolygon(DT.get_polygon_coordinates(vor, i, bbox)) for i in DT.each_polygon_index(vor)] )
end 
