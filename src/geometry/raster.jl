const __RASTER_TYPES = Union{SDMLayer, Matrix{<:Real}}


# ================================================================================
# TODO: 
# An alternative design decision here would be to not have a parameterized type,
# but to convert any matrix/other-rasterlike object to a SDMLayer.
# 
# A different alternative is enforce the SDT-useful fields here
    

"""
    Raster

A `Raster` stores gridded data. The internal representation of this data can be
an SDMLayer (from
[SimpleSDMLayers.jl](https://github.com/PoisotLab/SpeciesDistributionToolkit.jl))
or a Matrix.

`Raster` extends the RasterTrait type from GeoInterface.
"""
struct Raster{R <: __RASTER_TYPES}
    raster::R
end
Base.show(io::IO, r::Raster) = print(io, "Raster with dimensions $(size(r))")
Base.size(r::Raster) = size(r.raster)
Base.size(r::Raster, i::Integer) = size(r.raster, i)

Base.findall(r::Raster{<:SDMLayer}) = findall(r.raster.indices)
Base.findall(r::Raster{<:Matrix}) = findall(x->!ismissing(x) && !isnothing(x) && !isnan(x), r.raster)


Base.convert(::Type{Raster}, sdmlayer::SDMLayer) = Raster(sdmlayer)
Base.convert(::Type{Raster}, m::Matrix) = Raster(m)

Base.getindex(r::Raster, i::Integer) = getindex(r.raster, i)
Base.getindex(r::Raster, i::CartesianIndex) = getindex(r.raster, i)
Base.getindex(r::Raster, idx::Vector{<:CartesianIndex}) = map(i->getindex(r.raster, i), idx)
Base.getindex(r::Raster, bmat::BitMatrix) = [getindex(r.raster, i) for i in findall(bmat)]
Base.getindex(r::Raster, node::Node) = r.raster[node.coordinate...]
Base.getindex(r::Raster, bon::BiodiversityObservationNetwork) = [r.raster[node.coordinate...] for node in bon.nodes]

Base.eachindex(r::Raster) = eachindex(r.raster)

_get_cartesian_idx(r::Raster, node::Node) = CartesianIndex(SDT.SimpleSDMLayers.__get_grid_coordinate_by_latlon(r.raster, node.coordinate...))
_get_cartesian_idx(r::Raster, bon::BiodiversityObservationNetwork) = [_get_cartesian_idx(r, node) for node in bon]


Raster(sdmlayer::SDMLayer) = Raster{typeof(sdmlayer)}(sdmlayer)

datatype(::Raster{T}) where T = T.parameters[begin]
is_rasterizable(::T) where T = T <: __RASTER_TYPES


nonempty(r::Raster{<:SDMLayer}) = findall(r.raster.indices);
nonempty(r::Raster{<:Matrix}) = findall(x-> !isnothing(x) && !isnan(x) && !ismissing(x), r.raster)

_get_raster_extent(raster::Raster{<:SDMLayer}) = begin
    bbox = SDT.boundingbox(raster.raster)
    GI.Extent(X=(bbox.left, bbox.right), Y=(bbox.bottom, bbox.top))
end 
_get_raster_extent(::Raster{<:Matrix}) = GI.Extent(X=(0,1), Y=(0,1)) # rent on the unit square is out of control

_get_raster_crs(raster::Raster{<:SDMLayer}) = GI.crs(raster.raster)
_get_raster_crs(::Raster{<:Matrix}) = nothing


# SpeciesDistributionToolkit Overloads
SDT.eastings(r::Raster{<:Matrix}) = 0:(1/size(r)[2]):1
SDT.northings(r::Raster{<:Matrix}) = 0:(1/size(r)[1]):1
SDT.eastings(r::Raster{<:SDMLayer}) = eastings(r.raster)
SDT.northings(r::Raster{<:SDMLayer}) = northings(r.raster)

# GeoInterface overloads
GI.isgeometry(::Raster)::Bool = true
GI.geomtrait(::Raster)::DataType = GI.RasterTrait()
GI.israster(::Raster)::Bool = true
GI.trait(::Raster) = RasterTrait()
GI.extent(::RasterTrait, raster::Raster)::GI.Extents.Extent = _get_raster_extent(raster)
GI.crs(::RasterTrait, raster::Raster)::GeoFormatTypes.CoordinateReferenceSystem = _get_raster_crs(raster)