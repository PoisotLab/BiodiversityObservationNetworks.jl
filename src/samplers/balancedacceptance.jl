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

function _sample(sampler::BalancedAcceptance, geometry::T) where T<:Union{Polygon,Raster,Vector{<:Polygon}}
    _balanced_acceptance(sampler, geometry)
end 

function _sample(sampler::BalancedAcceptance, layers::RasterStack)
    _balanced_acceptance(sampler, first(layers))
end 

_get_easting_and_northing(::BalancedAcceptance, raster::Raster) = SDT.eastings(raster), SDT.northings(raster)
_get_easting_and_northing(sampler::BalancedAcceptance, layers::RasterStack) = _get_easting_and_northing(sampler, first(layers))
_get_easting_and_northing(sampler::BalancedAcceptance, polygon::Polygon) = begin 
    x, y = GI.extent(polygon)
    grid_size = sampler.grid_size
    Δx = (x[2]-x[1])/grid_size[1]
    Δy = (y[2]-y[1])/grid_size[2]

    Es = [x[1] + i*Δx for i in 1:grid_size[1]]
    Ns = [y[1] + i*Δy for i in 1:grid_size[2]]

    return Es, Ns
end 

# TODO: this is redundant and a similar thing is in BalancedAcceptance, unify
_check_candidate(Es, Ns, candidate, polygon::Polygon) = GeometryOps.contains(polygon, (Es[candidate[1]], Ns[candidate[2]]))
function _check_candidate(_, _, coord, raster::Raster)
    val = raster.raster[coord[2],coord[1]]
    !isnothing(val) && !ismissing(val) && !isnan(val)
end

function _balanced_acceptance(sampler, geometry) 
    num_nodes = sampler.number_of_nodes
    Es, Ns = _get_easting_and_northing(sampler, geometry)

    x_dim, y_dim = length(Es), length(Ns)

    seed = rand(Int.(1e0:1e7), 2)
    selected_points = Node[]
    ct = 0
    candct = 0
    while ct < num_nodes
        i, j = haltonvalue(seed[1] + candct, 2), haltonvalue(seed[2] + candct, 3)
        candct += 1
        candx, candy = convert.(Int, [ceil(x_dim * i), ceil(y_dim * j)])
        candidate = CartesianIndex(candx,candy)

        if _check_candidate(Es, Ns, candidate, geometry)
            push!(selected_points, Node((Es[candidate[1]], Ns[candidate[2]])))
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
