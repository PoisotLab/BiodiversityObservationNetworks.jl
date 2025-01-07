"""
    GeneralizedRandomTessellatedStratified

`GeneralizedRandomTessellatedStratified` is a type of [`BONSampler`](@ref) for
generating [`BiodiversityObservationNetwork`](@ref)s with spatial spreading.

*Arguments*:
- `number_of_nodes`: the number of sites to select
- `grid_size`: if being used on a polygon, the dimensions of the grid used to
  cover the extent. GRTS sampling uses discrete Cartesian indices


@Olsen

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
_valid_geometries(::GeneralizedRandomTessellatedStratified) = (Polygon, Raster, Vector{Polygon}, RasterStack)

"""
    _sample(sampler::GeneralizedRandomTessellatedStratified, geometry)

Internal dispatch for sampling using
[`GeneralizedRandomTessellatedStratified`](@ref) on a geometry.
"""
_sample(sampler::GeneralizedRandomTessellatedStratified, geometry::T) where T<:Union{Polygon,Raster,Vector{<:Polygon}} = _grts(sampler, geometry)
_sample(sampler::GeneralizedRandomTessellatedStratified, layers::RasterStack) = _grts(sampler, first(layers))


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
For `[Raster]`(@ref)s, this is computed based on the raster dimensions. 
"""
_get_address_length(sampler, ::Polygon) = Int(ceil(max(log(2, sampler.grid_size[1]), log(2, sampler.grid_size[2]))))
_get_address_length(sampler, raster::Raster) = Int(ceil(max(log(2, size(raster,1)), log(2, size(raster,2)))))

"""
    _get_easting_and_northing

Rn this is just copy/pasted from BAS. It should dispatch on a set of samplers,
with the assumption that the sampler always has a field called grid_size (this
also only matters for polygons)
"""
_get_easting_and_northing(::GeneralizedRandomTessellatedStratified, raster::Raster) = SDT.eastings(raster), SDT.northings(raster)
_get_easting_and_northing(sampler::GeneralizedRandomTessellatedStratified, layers::RasterStack) = _get_easting_and_northing(sampler, first(layers))
_get_easting_and_northing(sampler::GeneralizedRandomTessellatedStratified, polygon::Polygon) = begin
    x, y = GI.extent(polygon)
    grid_size = sampler.grid_size
    Δx = (x[2]-x[1])/grid_size[1]
    Δy = (y[2]-y[1])/grid_size[2]

    Es = [x[1] + i*Δx for i in 1:grid_size[1]]
    Ns = [y[1] + i*Δy for i in 1:grid_size[2]]
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

_get_cartesian_index_bounds(sampler, raster::Raster) = size(raster)
_get_cartesian_index_bounds(sampler, ::Polygon) = sampler.grid_size

_check_validity(geometry::Polygon, coord) = GO.contains(geometry, coord)

_check_validity(geometry::Raster, coord) = begin 
    val = geometry.raster[coord...]
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
            
            if _check_validity(geometry, coord)
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
#=
"""
    GeneralizedRandomTessellatedStratified

@Olsen
"""
@kwdef struct GeneralizedRandomTessellatedStratified <: BONSampler
    numsites = 50
    dims = (100, 100)
end

maxsites(grts::GeneralizedRandomTessellatedStratified) = prod(grts.dims)

function check_arguments(grts::GeneralizedRandomTessellatedStratified)
    check(TooManySites, grts)
    check(TooFewSites, grts)
    return
end

function _quadrant_fill!(mat)
    x, y = size(mat)
    a, b = x ÷ 2, y ÷ 2
    quad_ids = shuffle([1, 2, 3, 4])
    mat[begin:a, begin:b] .= quad_ids[1]
    mat[(a + 1):end, begin:b] .= quad_ids[2]
    mat[begin:a, (b + 1):end] .= quad_ids[3]
    mat[(a + 1):end, (b + 1):end] .= quad_ids[4]
    return mat
end

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
_generate!(coords::Vector{CartesianIndex}, grts::GeneralizedRandomTessellatedStratified)
"""
function _generate!(
    coords::Vector{CartesianIndex},
    grts::GeneralizedRandomTessellatedStratified,
)
    x, y = grts.dims # smallest multiple of 4 on each side
    num_address_grids = Int(ceil(max(log(2, x), log(2, y))))

    grid_sizes = reverse([2^i for i in 1:num_address_grids])

    address_grids =
        [zeros(Int, 2^num_address_grids, 2^num_address_grids) for _ in grid_sizes]

    map(
        i -> _quadrant_split!(address_grids[i], grid_sizes[i]),
        eachindex(address_grids),
    )

    code_numbers = sum([10^(i - 1) .* ag for (i, ag) in enumerate(address_grids)])
    sort_idx = sortperm([code_numbers[cidx] for cidx in eachindex(code_numbers)])

    return filter(
        idx -> idx[1] <= grts.dims[1] && idx[2] <= grts.dims[2],
        CartesianIndices(code_numbers)[sort_idx],
    )[1:(grts.numsites)]
end

=#