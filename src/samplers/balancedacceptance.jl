"""
    BalancedAcceptance

`BalancedAcceptance` is a type of [`BONSampler`](@ref) for generating
[`BiodiversityObservationNetwork`](@ref)s with spatial spreading.

First proposed in [Robertson2013BasBal](@cite). 

*Arguments*:
- `number_of_nodes`: the number of sites to select
- `grid_size`: if being used on a polygon, the dimensions of the grid used to
  cover the extent. Balanced Acceptance sampling uses discrete Cartesian indices
"""
Base.@kwdef struct BalancedAcceptance{I<:Integer} <: BONSampler
    number_of_nodes::I = 100
    grid_size::Tuple{I,I} = (250, 250)
end 
BalancedAcceptance(n::Integer; grid_size=(250, 250)) = BalancedAcceptance(n, grid_size)

_valid_geometries(::BalancedAcceptance) = (Polygon, Raster, Vector{Polygon}, RasterStack)

function _sample(sampler::BalancedAcceptance, geometry::T) where T<:Union{Polygon,SDMLayer,Vector{<:Polygon}}
    _balanced_acceptance(sampler, geometry)
end 

function _sample(sampler::BalancedAcceptance, layers::Vector{<:SDMLayer})
    _balanced_acceptance(sampler, first(layers))
end 

_get_easting_and_northing(::BalancedAcceptance, raster::SDMLayer) = SDT.eastings(raster), SDT.northings(raster)
_get_easting_and_northing(sampler::BalancedAcceptance, layers::Vector{<:SDMLayer}) = _get_easting_and_northing(sampler, first(layers))
_get_easting_and_northing(sampler::BalancedAcceptance, polygon::Polygon) = begin 
    x, y = GI.extent(polygon)
    grid_size = sampler.grid_size

    # these are flipped to match the behavior of Rasters
    easting_ticks, northing_ticks = grid_size[2], grid_size[1]

    Δx = (x[2]-x[1])/easting_ticks
    Δy = (y[2]-y[1])/northing_ticks

    Es = [x[1] + i*Δx for i in 1:easting_ticks]
    Ns = [y[1] + i*Δy for i in 1:northing_ticks]
    return Es, Ns
end 

# TODO: this is redundant and a similar thing is in BalancedAcceptance, unify
_check_candidate(Es, Ns, candidate, polygon::Polygon) = GeometryOps.contains(polygon, (Es[candidate[2]], Ns[candidate[1]]))
function _check_candidate(_, _, coord, raster::SDMLayer)
    (coord[1] > size(raster, 1) || coord[2] > size(raster, 2)) && return false
    val = raster[coord[1],coord[2]]
    !isnothing(val) && !ismissing(val) && !isnan(val)
end

function _balanced_acceptance(sampler, geometry) 
    num_nodes = sampler.number_of_nodes

    Es, Ns = _get_easting_and_northing(sampler, geometry)
    num_eastings, num_northings = length(Es), length(Ns)

    seed = rand(Int.(1e0:1e7), 2)
    selected_points = Node[]
    ct = 0
    candct = 0
    while ct < num_nodes
        i, j = haltonvalue(seed[1] + candct, 2), haltonvalue(seed[2] + candct, 3)
        candct += 1

        # northings are the first dim, eastings are the second
        candidate = CartesianIndex(convert.(Int, [ceil(num_northings * i), ceil(num_eastings * j)])...)
        if _check_candidate(Es, Ns, candidate, geometry)
            push!(selected_points, Node((Es[candidate[2]], Ns[candidate[1]])))
            ct += 1
        end
    end
    return BiodiversityObservationNetwork(selected_points)
end 

function _sample(sampler::BalancedAcceptance, domain::Vector{<:Polygon})
    @info "You passed a Vector of Polygons."
    @info "Note by default BalancedAcceptance applies to each Polygon separately"
    @info "To use BalancedAcceptance on the whole extent, merge the polygons."

    vcat([_sample(sampler, p) for p in domain])
end 

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use BalancedAcceptance with default arguments on a Raster" begin
    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(zeros(50, 100))
    bas = BalancedAcceptance()
    bon = sample(bas, raster)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == bas.number_of_nodes
end


@testitem "We can use BalancedAcceptance with default arguments on a Polygon" begin
    polygon = openstreetmap("COL")
    bas = BalancedAcceptance()
    bon = sample(bas, polygon)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == bas.number_of_nodes
end