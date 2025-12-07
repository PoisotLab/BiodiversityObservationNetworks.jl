"""
    to_domain(x; kwargs...)

Coerce an input object into a valid internal domain representation (e.g., `RasterDomain`, 
`RasterStack`).

# Arguments
- `x`: The input object to convert. Can be a Matrix, SDMLayer, Vector of layers, 
  Polygon, or an existing domain/network.
- `kwargs...`: Helper arguments for specific conversions (see below).

# Methods

- **1. Matrices and SDMLayers**: Wraps the input in a RasterDomain. For SDMLayer inputs, the existing valid indices (non-NaN/masked) are preserved in the pool.
    - `to_domain(mat::AbstractMatrix)` -> `RasterDomain`
    - `to_domain(layer::SDMLayer)` -> `RasterDomain`

- **2. Collections of layers (Stacks)**: Converts a vector of inputs into a RasterStack. The sampling pool is aggregated as the intersection of valid pixels across all layers.
    - `to_domain(layers::Vector{<:SDMLayer})` -> `RasterStack`
    - `to_domain(mats::Vector{<:AbstractMatrix})` -> `RasterStack`
    
- **3. Polygons (Rasterization)**: Converts a vector geometry (Polygon) into a binary RasterDomain by rasterizing it.
    - `to_domain(poly::SimpleSDMPolygons.AbstractGeometry; grid_size=(100,100))` -> `RasterDomain`
    - `to_domain(poly::PolygonDomain; grid_size=(100,100))` -> `RasterDomain`

- **4. Identity**: Returns the input unchanged if it is already a valid domain or network type.
    - `to_domain(rd::RasterDomain)` -> `RasterDomain`
    - `to_domain(rs::RasterStack)` -> `RasterStack`
    - `to_domain(bon::BiodiversityObservationNetwork)` -> `BiodiversityObservationNetwork`

"""
to_domain(layer::SDMLayer; kwargs...) = RasterDomain(layer, layer.indices)
to_domain(mat::AbstractMatrix; kwargs...) = RasterDomain(Matrix(mat))

# Collections of layers
to_domain(layers::Vector{<:SDMLayer}; kwargs...) = RasterStack(layers)
to_domain(rds::Vector{<:RasterDomain}; kwargs...) = RasterStack(rds)
to_domain(mat::Vector{<:AbstractMatrix}; kwargs...) = RasterStack(mat)

# Polygons
to_domain(poly::SimpleSDMPolygons.AbstractGeometry; kwargs...) = to_domain(PolygonDomain(poly); kwargs...)
function to_domain(poly::PolygonDomain; grid_size=(100,100))
    # Create blank SDMLayer
    x,y = extent(poly)
    layer = RasterDomain(SDMLayer(ones(Bool, grid_size), x=x, y=y))

    # Mask it by polygon
    mask!(layer, poly)
    return layer
end

# Identity 
to_domain(rd::RasterDomain; kwargs...) = rd
to_domain(rs::RasterStack, kwargs...) = rs
to_domain(bon::BiodiversityObservationNetwork; kwargs...) = bon

# ----------------------------------------------------------------
# Conversion to polygons. 
# Not sure if this is needed anymore
# ----------------------------------------------------------------
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
