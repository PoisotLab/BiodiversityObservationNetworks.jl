struct RasterStack{R<:Raster}
    stack::Vector{R}
end
Base.show(io::IO, ls::RasterStack) = print(io, "RasterStack with $(length(ls.stack)) layers")
Base.getindex(ls::RasterStack, i::Integer)= ls.stack[i]
Base.length(layers::RasterStack) = length(layers.stack)
Base.iterate(layers::RasterStack, i) = iterate(layers.stack, i)
Base.iterate(layers::RasterStack) = iterate(layers.stack)

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


features(layers::RasterStack) = findall(layers[1].raster.indices), hcat([layer.raster.grid[layer.raster.indices] for layer in layers]...)'
