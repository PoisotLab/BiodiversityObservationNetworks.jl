# ------------------------------------------
#  Domain normalization 
#
# ------------------------------------------
to_domain(layer::SDMLayer; kwargs...) = RasterDomain(layer, layer.indices)
to_domain(rd::RasterDomain; kwargs...) = rd
to_domain(mat::AbstractMatrix; kwargs...) = RasterDomain(Matrix(mat))


to_domain(layers::Vector{<:SDMLayer}; kwargs...) = RasterStack(layers)
to_domain(rds::Vector{<:RasterDomain}; kwargs...) = RasterStack(rds)
to_domain(mat::Vector{<:AbstractMatrix}; kwargs...) = RasterStack(mat)

to_domain(rs::RasterStack, kwargs...) = rs


to_domain(poly::SimpleSDMPolygons.AbstractGeometry; kwargs...) = to_domain(PolygonDomain(poly); kwargs...)

function to_domain(poly::PolygonDomain; grid_size=(100,100))
    # Create blank SDMLayer
    x,y = extent(poly)
    layer = RasterDomain(SDMLayer(ones(Bool, grid_size), x=x, y=y))

    # Mask it by polygon
    mask!(layer, poly)
    return layer
end


to_domain(bon::BiodiversityObservationNetwork; kwargs...) = bon


function to_polygon(domain::RasterDomain)
    x, y = extent(domain)
    x, y = Float32.(x), Float32.(y)
    bbox = (left=x[1], right=x[2], bottom=y[1], top=y[2])
    poly = SimpleSDMPolygons._get_polygon_from_bbox(bbox)
    return PolygonDomain(poly)
end

to_polygon(mat::Matrix) = to_polygon(to_domain(mat))
to_polygon(layer::SDMLayer) = to_polygon(to_domain(layer))
to_polygon(poly::SimpleSDMPolygons.AbstractGeometry) = PolygonDomain(poly)
to_polygon(poly::PolygonDomain) = poly
