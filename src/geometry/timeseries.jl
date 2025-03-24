
function _common_mask!(sdmlayers::Vector{S}) where S<:SDMLayer
    mask_grid = reduce(.&, [l.indices for l in sdmlayers])
    for layer in sdmlayers
        layer.indices .= mask_grid
    end 
end

function _common_mask!(timeseries::Vector{<:Vector{S}}) where S<:SDMLayer
    idx = vcat([[r.indices for r in st] for st in timeseries]...)
    mask_grid = reduce(.&, idx)
    for stack in timeseries
        for r in stack 
            r.indices .= mask_grid
        end 
    end 
end 

_baseline_worldclim_timespan() = Year(1970)=>Year(2000)
_worldclim_timespans() = [Year(y)=>Year(y+19) for y in [2021, 2041, 2061, 2081]]

function _get_bioclim_baseline(
    dataprovider,
    mask,
    extent
)
    ex = isnothing(extent) ? (isnothing(mask) ? () : SDT.boundingbox(mask)) : extent
    bioclim = [SDT.SDMLayer(dataprovider; layer="BIO$i", ex...) for i in 1:19]
    !isnothing(mask) && SDT.mask!(bioclim, mask)
    return bioclim
end

function _get_single_bioclim_future(
    dataprovider,
    future,
    timespan, 
    mask,
    extent,
    idx
)
    function _get_layer(number)
        ex = isnothing(extent) ? (isnothing(mask) ? () : SDT.boundingbox(mask)) : extent
        layer = SDMLayer(dataprovider, future; layer="BIO$number", timespan=timespan, ex...)
        layer.indices .= idx
        layer.indices[findall(isnan, layer.grid)] .= 0  
        return layer
    end 
    
    layers = _get_layer.(1:19)
    !isnothing(mask) && SDT.mask!(layers, mask)
    return layers
end


"""
    bioclim_futures

A convenience method for getting the BioClim variable time-series under different SSPs.
"""
function bioclim_futures(
    ; 
    ssp = SDT.SSP245,
    mask = nothing,
    extent = nothing, 
    earth_system_model = SDT.ACCESS_CM2
)

    timespans = _worldclim_timespans()

    dataprovider = RasterData(WorldClim2, BioClim)
    future = Projection(ssp, earth_system_model)

    baseline = _get_bioclim_baseline(dataprovider, mask, extent)
    idx = baseline[1].indices
    stacks = vcat([baseline], [_get_single_bioclim_future(dataprovider, future, ts, mask, extent, idx) for ts in timespans])
    _common_mask!(stacks)
    return stacks
end


