"""
    RasterDomain{T, P} <: AbstractDomain

A wrapper around a raster-like object (Matrix or SDMLayer) that maintains a "pool"
of valid sampling indices.

# Fields
- `data::T`: The underlying raster data (e.g., `Matrix` or `SDMLayer`).
- `pool::P`: A collection of valid indices that can be sampled.
"""
struct RasterDomain{T}
    data::T
    pool
end

RasterDomain(data) = RasterDomain(data, ones(Bool, size(data)))

# Base overloads for RasterDomain
Base.size(raster::RasterDomain) = size(raster.data)
Base.size(raster::RasterDomain, i) = size(raster.data, i)
Base.length(raster::RasterDomain) = length(raster.data)

Base.iterate(raster::RasterDomain) = iterate(raster.data)
Base.iterate(raster::RasterDomain, i) = iterate(raster.data, i)

Base.getindex(raster::RasterDomain, i) = getindex(raster.data[i])
Base.getindex(raster::RasterDomain, i, j) = getindex(raster.data, i,j)
Base.getindex(raster::RasterDomain, ci::CartesianIndex) = raster.data[ci]
Base.getindex(raster::RasterDomain, i::Tuple) = raster[i...]
Base.setindex!(raster::RasterDomain, v, i) = setindex!(raster.data, v, i)

Base.sum(raster::RasterDomain) = sum(raster.data)

"""
    getpool(domain)

Return the collection of valid (unmasked) indices for the given domain.
"""
getpool(raster::RasterDomain) = findall(raster.pool)


"""
    getfeatures(domain)

Return a matrix of auxiliary variables (features) for valid pixels in the domain.
Rows are features, columns are pixels.
"""
function getfeatures(rd::RasterDomain)
    pool = getpool(rd)
    return hcat([rd[i] for i in pool]...)
end

"""
    getcoordinates(domain)

Retrieve the spatial coordinates of all valid (unmasked) sampling locations in the domain.

# Arguments
- `domain`: A `RasterDomain` or a `RasterStack`.

# Returns
- A `2 x N` Matrix{Float32}, where `N` is the number of valid locations in the pool, 
where each column is a coordinate in the valid pool of locations.

# Description
Returns a matrix of coordinates corresponding to the valid sampling pool. 
- For **Matrix-backed domains**, coordinates are the integer column (x) and row (y) indices.
- For **SDMLayer-backed domains**, coordinates are the projected spatial coordinates 
  (e.g., Longitude/Latitude) derived from the layer's geotransform.
- For **RasterStacks**, coordinates are derived from the first layer in the stack 
  (assuming all layers share the same grid).
"""
function getcoordinates(raster::RasterDomain{<:Matrix})
    return Float32.(hcat([[x[1], x[2]] for x in getpool(raster)]...))
end

function getcoordinates(raster::RasterDomain{<:SDMLayer})
    Es, Ns = SpeciesDistributionToolkit.eastings(raster.data), SpeciesDistributionToolkit.northings(raster.data)
    return hcat([[j,i] for i in Ns, j in Es][raster.pool]...)
end


"""
    extent
"""
extent(rd::RasterDomain{<:Matrix}) = Extent(X=(1,size(rd,1)), Y=(1, size(rd, 2)))
extent(rd::RasterDomain{<:SDMLayer}) = Extent(X=rd.data.x, Y= rd.data.y)
extent(layer::SDMLayer) = Extent(X=layer.x, Y= rd.y)

"""
    crs
"""
crs(::RasterDomain{<:Matrix}) = nothing
crs(rd::RasterDomain{<:SDMLayer}) = crs(rd.data)
crs(layer::SDMLayer) = layer.crs

"""
    convert_node
"""
convert_node(::RasterDomain{<:Matrix}, cidx) = cidx
convert_node(raster::RasterDomain{<:SDMLayer}, cidx) = eastings(raster.data)[cidx[2]], northings(raster.data)[cidx[1]] 

"""
    convert_nodes
"""
convert_nodes(raster::RasterDomain, nodes) = [convert_node(raster, n) for n in nodes]


"""
    rescale_node
"""
function rescale_node(domain::RasterDomain{<:Matrix}, x::Real, y::Real) 
    x_scaled, y_scaled = Int.(ceil.(size(domain) .* [x, y]))
    return x_scaled, y_scaled
end    

function rescale_node(domain::RasterDomain{<:SDMLayer}, x_offset::Real, y_offset::Real) 
    (xmin, xmax), (ymin, ymax) = extent(domain)
    x_scaled = xmin + (xmax - xmin)*x_offset
    y_scaled = ymin + (ymax - ymin)*y_offset
    return x_scaled, y_scaled
end     

"""
    ismasked
"""
function ismasked(rd::RasterDomain{<:SDMLayer}, i::Real, j::Real) 
    x, y = SimpleSDMLayers.__get_grid_coordinate_by_latlon(rd.data, i, j)
    return !rd.pool[x,y]
end 

ismasked(rd::RasterDomain{<:SDMLayer}, i::Integer, j::Integer) = ismasked(rd, CartesianIndex(i,j))
ismasked(rd::RasterDomain, i::Integer, j::Integer) = ismasked(rd, CartesianIndex(i,j))

ismasked(rd::RasterDomain, ci::CartesianIndex) = !rd.pool[ci]

