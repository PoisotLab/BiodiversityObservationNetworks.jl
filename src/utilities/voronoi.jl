
function voronoi(bon::BiodiversityObservationNetwork, geom)
    tri = DT.triangulate([n.coordinate for n in bon])
    vor = DT.voronoi(tri)

    (xmin, xmax), (ymin, ymax) = GI.extent(geom) 
    bbox = (xmin, xmax, ymin, ymax) 


    # clipped_coords = Vector{Vector{NTuple{2, Float64}}}(undef, DT.num_polygons(vor))
    #for i in DT.each_polygon_index(vor)
    #    clipped_coords[i] = DT.get_polygon_coordinates(vor, i, bbox)
    #end
    return convert.(Polygon, [AGDAL.createpolygon(DT.get_polygon_coordinates(vor, i, bbox))
    for i in DT.each_polygon_index(vor)] )


    agdal_poly = AGDAL.createmultipolygon(polys) 
    return convert(Polygon, agdal_poly)
end 