"""
    Grid

`Grid` is a type of [`BONSampler`](@ref) for generating
[`BiodiversityObservationNetwork`](@ref)s with [`Node`](@ref)s structured as a
systematic grid over the study area.

*Arguments*:
- `longitude_spacing`
- `latitude_spacing`
"""
Base.@kwdef struct Grid{F<:Real} <: BONSampler
    longitude_spacing::F = 1. # in wgs84 coordinates
    latitude_spacing::F  = 1.
    # there should probably be padding here. 
end 

_valid_geometries(::Grid) = (Polygon, SDMLayer, Vector{<:SDMLayer})

_sample(sampler::Grid, layers::Vector{<:SDMLayer}) = _sample(sampler, first(layers))
function _sample(sampler::Grid, raster::SDMLayer)
    (xm, xM), (ym, yM) = GI.extent(raster)
    x_step, y_step = sampler.longitude_spacing, sampler.latitude_spacing
    BiodiversityObservationNetwork(vec([Node((i,j)) for i in xm:x_step:xM, j in ym:y_step:yM]))
end

function _sample(sampler::Grid, polygon::Polygon)
    (xm, xM), (ym, yM) = GI.extent(polygon)
    x_step, y_step = sampler.longitude_spacing, sampler.latitude_spacing
    BiodiversityObservationNetwork(vec([Node((i,j)) for i in xm:x_step:xM, j in ym:y_step:yM if GO.contains(polygon, (i,j))]))
end 


# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use a Grid with default constructor on a Raster" begin
    gs = Grid()
    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(zeros(50,30))
    bon = sample(gs, raster)
    @test bon isa BiodiversityObservationNetwork
end


@testitem "We can use a Grid with default constructor on a Polygon" begin
    gs = Grid()
    poly = openstreetmap("COL")
    bon = sample(gs, poly)
    @test bon isa BiodiversityObservationNetwork
end