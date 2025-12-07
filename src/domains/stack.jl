"""
    RasterStack{T} <: AbstractDomain

A collection of aligned `RasterDomain`s, used when sampling requires multivariate features.

# Fields
- `rasters::Vector{<:RasterDomain{<:T}}`: The stack of layers.
- `pool`: The intersection of valid pools from all layers.
"""
struct RasterStack{T}
    rasters::Vector{<:RasterDomain{<:T}}
    pool
end

# Primary constructor from a vector of `RasterDomain`s
function RasterStack(layers::Vector{<:T}) where T<:Union{RasterDomain,SDMLayer,AbstractMatrix}
    domains = to_domain.(layers)

    # assert same extent, crs, size for each layer
    allequal(crs.(domains)) || throw(ArgumentError("Not all layers have the same CRS"))
    allequal(size.(domains)) || throw(ArgumentError("Not all layers have the same size"))
    #allequal(extent.(domains)) || throw(ArgumentError("Not all layers have the same extent"))

    total_pool = _aggregate_pool(domains)
    RasterStack(domains, total_pool)
end


# Base overloads 
Base.show(io::IO, rs::RasterStack{T}) where T = print(io, "RasterStack{$T} with $(length(rs.rasters)) layers")
Base.size(rs::RasterStack) = size(first(rs.rasters))
Base.size(rs::RasterStack, i) = size(first(rs.rasters), i)
Base.length(rs::RasterStack) = length(first(rs.rasters))

Base.iterate(rs::RasterStack) = iterate(rs.rasters)
Base.iterate(rs::RasterStack, i) = iterate(rs.rasters, i)

Base.getindex(rs::RasterStack, i) = [r[i] for r in rs.rasters]
Base.getindex(rs::RasterStack, i, j) = [r[i,j] for r in rs.rasters]

getpool(rs::RasterStack) = findall(first(rs.rasters).pool)

function getfeatures(rs::RasterStack)
    pool = getpool(rs)
    return hcat([rs[i] for i in pool]...)
end


function getcoordinates(rs::RasterStack)
    getcoordinates(rs[1])
end

"""
    _aggregate_pool(domains)

Compute the intersection of valid sampling pools across multiple domains.

# Arguments
- `domains`: An iterable collection of domain objects (e.g., a vector of `RasterDomain`). 
  Each element must have a `pool` field containing a boolean mask.

# Returns
- A `BitArray` representing the intersection of all pools. A cell 
  is considered valid (`true`) only if it is valid in *all* provided domains.

# Description
This internal utility is primarily used when constructing a `RasterStack`. It ensures 
that sampling only occurs in locations that are valid across every layer in the stack 
(e.g., avoiding `NaN`s, missing data, or masked areas present in any single layer).
"""
function _aggregate_pool(domains)
    pools = [d.pool for d in domains]
    total_pool = reduce((A, B) -> A .& B, pools) 
    return total_pool
end 


extent(rs::RasterStack) = extent(first(rs.rasters))
crs(rs::RasterStack) = crs(first(rs.rasters))

convert_node(rs::RasterStack, x...) = convert_node(first(rs.rasters), x...)
convert_nodes(rs::RasterStack, nodes) = convert_nodes(first(rs.rasters), nodes)
rescale_node(rs::RasterStack, x...) = rescale_node(first(rs.rasters), x...)

ismasked(rs::RasterStack, x...) = ismasked(first(rs.rasters), x...)

