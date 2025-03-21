struct RasterTimeseries{T<:Union{RasterStack,Raster},D<:Union{Pair{<:DatePeriod,<:DatePeriod},DateTime}}
    timeseries::Vector{T}
    dates::Vector{D}
    function RasterTimeseries(timeseries::Vector{<:I}, dates::Vector{<:J}) where {I,J}
        _common_mask!(timeseries)
        new{I,J}(timeseries, dates)
    end
end

_common_mask!(timeseries::Vector{<:Raster}) = _common_mask!(map(r->r.raster, timeseries))

function _common_mask!(timeseries::Vector{<:RasterStack}) 
    idx = vcat([[r.raster.indices for r in st] for st in timeseries]...)
    mask_grid = reduce(.&, idx)
    for stack in timeseries
        for r in stack 
            r.raster.indices .= mask_grid
        end 
    end 
end 


Base.show(io::IO, rt::RasterTimeseries) = print(io, "RasterTimeseries for $(length(first(rt.timeseries))) layers across $(length(rt.timeseries)) timesteps")
Base.getindex(rt::RasterTimeseries, i) = rt.timeseries[i]
Base.iterate(rt::RasterTimeseries) = iterate(rt.timeseries)
Base.iterate(rt::RasterTimeseries, i) = iterate(rt.timeseries, i)
Base.length(rt::RasterTimeseries) = length(rt.timeseries)

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
    return RasterStack(bioclim)
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
    RasterStack(Raster.(layers))
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
    idx = baseline[1].raster.indices
    stacks = vcat(baseline, [_get_single_bioclim_future(dataprovider, future, ts, mask, extent, idx) for ts in timespans])

    return RasterTimeseries(stacks, vcat(_baseline_worldclim_timespan(), timespans...))
end


