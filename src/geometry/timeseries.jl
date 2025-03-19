struct Timespan
    start::DateTime
    finish::DateTime
end


struct RasterTimeseries{T<:Union{RasterStack,Raster},D<:Union{Pair{<:DatePeriod,<:DatePeriod},DateTime}}
    timeseries::Vector{T}
    dates::Vector{D}
end
Base.show(io::IO, rt::RasterTimeseries) = print(io, "RasterTimeseries for $(length(first(rt.timeseries))) layers across $(length(rt.timeseries)) timesteps")


_baseline_worldclim_timespan() = Year(1970)=>Year(2000)
_worldclim_timespans() = [Year(y)=>Year(y+19) for y in [2021, 2041, 2061, 2081]]

function _get_bioclim_baseline(
    dataprovider,
    mask,
    extent
)
    ex = isnothing(extent) ? () : extent
    bioclim = [SDT.SDMLayer(dataprovider; layer="BIO$i", ex...) for i in 1:19]
    !isnothing(mask) && SDT.mask!(layers, mask)
    return RasterStack(bioclim)
end

function _get_single_bioclim_future(
    dataprovider,
    future,
    timespan, 
    mask,
    extent,
)
    function _get_layer(number)
        ex = isnothing(extent) ? () : extent
        layer = SDMLayer(dataprovider, future, layer="BIO$number", timespan=timespan, ex...)
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

    stacks = vcat(_get_bioclim_baseline(dataprovider, mask, extent), [_get_single_bioclim_future(dataprovider, future, ts, mask, extent) for ts in timespans])

    return RasterTimeseries(stacks, vcat(_baseline_worldclim_timespan(), timespans...))
end


