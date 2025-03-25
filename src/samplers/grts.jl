"""
    GeneralizedRandomTessellatedStratified

`GeneralizedRandomTessellatedStratified` is a type of [`BONSampler`](@ref) for
generating [`BiodiversityObservationNetwork`](@ref)s with spatial spreading.

GRTS was initially proposed in [Stevens2004SpaBal](@cite).

*Arguments*:
- `number_of_nodes`: the number of sites to select
- `grid_size`: if being used on a polygon, the dimensions of the grid used to
  cover the extent. GRTS sampling uses discrete Cartesian indices


GRTS represents each cell of a rasterized version of the sampling domain
using an address, where the address of each cell is represented as a `D`-digit
base-4 number. 

The value of `D` depends on the size of the raster. GRTS works by recursively
splitting the rasterized domain into quadrants, and those quadrants into further
quadrants, until a single pixel is reached. At each step, each quadrant is
randomly labeled with a number 1, 2, 3, or 4 (without replacement, so all values are
used). 

The addresses are then sorted numerically, and the `number_of_nodes` smallest
values are chosen.

"""
Base.@kwdef struct GeneralizedRandomTessellatedStratified{I<:Integer} <: BONSampler
    number_of_nodes::I = 100
    grid_size::Tuple{I,I} = (250, 250)
end 
GeneralizedRandomTessellatedStratified(n::Integer; grid_size=(250,250)) = GeneralizedRandomTessellatedStratified(n, grid_size)


_valid_geometries(::GeneralizedRandomTessellatedStratified) = (Polygon, SDMLayer, Vector{Polygon}, Vector{<:SDMLayer})

"""
    _sample(sampler::GeneralizedRandomTessellatedStratified, geometry)

Internal dispatch for sampling using
[`GeneralizedRandomTessellatedStratified`](@ref) on a geometry.
"""
_sample(sampler::GeneralizedRandomTessellatedStratified, geometry::T) where T<:Union{Polygon,SDMLayer,Vector{<:Polygon}} = _grts(sampler, geometry)
_sample(sampler::GeneralizedRandomTessellatedStratified, layers::Vector{<:SDMLayer}) = _grts(sampler, first(layers))


"""
    _quadrant_fill!(mat)

Takes a matrix `mat` and splits it into quadrants randomly labeled one through four. 
"""
function _quadrant_fill!(mat)
    x, y = size(mat)
    a, b = x ÷ 2, y ÷ 2
    quad_ids = Random.shuffle([1, 2, 3, 4])
    mat[begin:a, begin:b] .= quad_ids[1]
    mat[(a + 1):end, begin:b] .= quad_ids[2]
    mat[begin:a, (b + 1):end] .= quad_ids[3]
    mat[(a + 1):end, (b + 1):end] .= quad_ids[4]
    return mat
end


"""
    _quadrant_split!(mat, grid_size)

Splits a matrix `mat` into nested quadrants, where the side-length of a
submatrix to be split into quadrants is given by `grid_size`.  
"""
function _quadrant_split!(mat, grid_size)
    x, y = size(mat)
    num_x_grids, num_y_grids = x ÷ grid_size, y ÷ grid_size

    for i in 0:(num_x_grids - 1), j in 0:(num_y_grids - 1)
        a, b = i * grid_size, j * grid_size
        bounds = (a + 1, b + 1), (a + grid_size, b + grid_size)
        mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]] .=
            _quadrant_fill!(mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]])
    end

    return mat
end

"""
    _get_address_length

Returns the number of digits in a [`GeneralizedRandomTessellatedStratified`](@ref)
address for a specific geometry. For [`Polygon`](@ref)s, this is computed based
on the `grid_size` field of the `GeneralizedRandomTessellatedStratified` object.
For SDMLayers, this is computed based on the raster dimensions. 
"""
_get_address_length(sampler, ::Polygon) = Int(ceil(max(log(2, sampler.grid_size[1]), log(2, sampler.grid_size[2]))))
_get_address_length(sampler, raster::SDMLayer) = Int(ceil(max(log(2, size(raster,1)), log(2, size(raster,2)))))

"""
    _get_easting_and_northing

Rn this is just copy/pasted from BAS. 

TODO: everything is fucked because Raster dims are represented as (northings/eastings)

It should dispatch on a set of samplers,
with the assumption that the sampler always has a field called grid_size (this
also only matters for polygons)
"""
_get_easting_and_northing(::GeneralizedRandomTessellatedStratified, raster::SDMLayer) = SDT.eastings(raster), SDT.northings(raster)
_get_easting_and_northing(sampler::GeneralizedRandomTessellatedStratified, layers::Vector{<:SDMLayer}) = _get_easting_and_northing(sampler, first(layers))
_get_easting_and_northing(sampler::GeneralizedRandomTessellatedStratified, polygon::Polygon) = begin
    #=
    x, y = GI.extent(polygon)
    grid_size = sampler.grid_size
    Δx = (x[2]-x[1])/grid_size[1]
    Δy = (y[2]-y[1])/grid_size[2]

    Es = [x[1] + i*Δx for i in 1:grid_size[1]]
    Ns = [y[1] + i*Δy for i in 1:grid_size[2]]
    return Es, Ns=#
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

"""
    _get_addresses()
"""
function _get_addresses(sampler, geometry)
    address_length = _get_address_length(sampler, geometry)
    grid_sizes = reverse([2^i for i in 1:address_length])
    grids = [zeros(Int, 2^address_length, 2^address_length) for _ in grid_sizes]

    map(
        i -> _quadrant_split!(grids[i], grid_sizes[i]),
        eachindex(grids),
    )
    addresses = sum([10^(i - 1) .* ag for (i, ag) in enumerate(grids)])
    return addresses
end

_get_cartesian_index_bounds(sampler, raster::SDMLayer) = size(raster)
_get_cartesian_index_bounds(sampler, ::Polygon) = sampler.grid_size

_check_validity(geometry::Polygon, coord) = GO.contains(geometry, coord)

_check_validity(geometry::SDMLayer, coord) = begin 
    val = geometry[coord[1], coord[2]]
    !isnothing(val) && !ismissing(val) && !isnan(val)
end 
function _pick_nodes(sampler, geometry, addresses)
    Es, Ns = _get_easting_and_northing(sampler, geometry)
    xbound, ybound = _get_cartesian_index_bounds(sampler, geometry)
    sort_idx = sortperm([addresses[cidx] for cidx in eachindex(addresses)])
    cart_idxs = CartesianIndices(addresses)[sort_idx]
    
    num_selected = 0
    cursor = 1

    selected_nodes = Node[]
    while num_selected < sampler.number_of_nodes && cursor <= length(sort_idx)
        candidate = cart_idxs[cursor]
        if candidate[1] <= xbound && candidate[2] <= ybound
            coord = (Es[candidate[2]], Ns[candidate[1]])
            
            # hacky, but it works
            c = geometry isa Polygon ? coord : candidate
            if _check_validity(geometry, c)
                push!(selected_nodes, Node(coord))
                num_selected += 1
            end
        end 
        cursor += 1
    end 
    return selected_nodes 
end

function _grts(sampler, geometry)
    addresses = _get_addresses(sampler, geometry)
    selected_nodes = _pick_nodes(sampler, geometry, addresses)
    return BiodiversityObservationNetwork(selected_nodes)
end 


# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use GRTS with default arguments on a Raster" begin
    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(zeros(50, 100))
    grts = GeneralizedRandomTessellatedStratified()
    bon = sample(grts, raster)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == grts.number_of_nodes
end


@testitem "We can use GRTS with default arguments on a Polygon" begin
    poly = gadm("COL")
    grts = GeneralizedRandomTessellatedStratified()
    bon = sample(grts, poly)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == grts.number_of_nodes
end