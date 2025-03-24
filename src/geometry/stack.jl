"""
    RasterStack

A `RasterStack` is a set of [`Raster`](@ref)s, which all have the same
resolution and extent, for storing multiple raster-like environmental variables 
"""
struct RasterStack{R<:Raster}
    stack::Vector{R}
end
Base.show(io::IO, ls::RasterStack) = print(io, "RasterStack with $(length(ls.stack)) layers")
Base.getindex(ls::RasterStack, i::Vector{<:Integer})= RasterStack(ls.stack[i])
Base.getindex(ls::RasterStack, i::Integer)= ls.stack[i]
Base.getindex(ls::RasterStack, i::CartesianIndex) = [getindex(r.raster, i) for r in ls]
Base.getindex(ls::RasterStack, idxs::Vector{<:CartesianIndex}) = Matrix(hcat([getindex(r, idxs) for r in ls]...)')

Base.getindex(rs::RasterStack, bon::BiodiversityObservationNetwork) = hcat([[r[node] for r in rs.stack] for node in bon.nodes]...)
Base.getindex(rs::RasterStack, node::Node) = [r.raster[node.coordinate...] for r in rs.stack]

Base.first(layers::RasterStack) = layers[1]
Base.size(layers::RasterStack) = size(first(layers))
Base.length(layers::RasterStack) = length(layers.stack)
Base.iterate(layers::RasterStack, i) = iterate(layers.stack, i)
Base.iterate(layers::RasterStack) = iterate(layers.stack)
Base.eachindex(layers::RasterStack) = Base.oneto(length(layers))


_get_cartesian_idx(layers::RasterStack, node::Node) = _get_cartesian_idx(first(layers), node)
_get_cartesian_idx(layers::RasterStack, bon::BiodiversityObservationNetwork) = _get_cartesian_idx(first(layers), bon)


function _common_mask!(sdmlayers::Vector{S}) where S<:SDMLayer
    mask_grid = reduce(.&, [l.indices for l in sdmlayers])
    for layer in sdmlayers
        layer.indices .= mask_grid
    end 
end 

function RasterStack(sdmlayers::Vector{S}) where S<:SDMLayer
    _common_mask!(sdmlayers)
    RasterStack(Raster.(sdmlayers))
end

SDT.eastings(layers::RasterStack) = eastings(first(layers))
SDT.northings(layers::RasterStack) = northings(first(layers))
nonempty(r::RasterStack) = ∩(nonempty.(r)...)

features(layers::Vector{<:SDT.SDMLayer}) = hcat([[l[x] for l in layers] for x in eachindex(first(layers))]...)
features(layers::RasterStack) = findall(layers[1].raster.indices), Matrix(hcat([layer.raster.grid[layer.raster.indices] for layer in layers]...)')
