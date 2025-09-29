"""
    GeneralizedRandomTesselated

`GeneralizedRandomTesselated` is a type of [`BONSampler`](@ref) for
generating [`BiodiversityObservationNetwork`](@ref)s with spatial spreading.

GRTS was initially proposed in [Stevens2004SpaBal](@cite).

*Arguments*:
- `num_nodes`: the number of sites to select
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

The addresses are then sorted numerically, and the `num_nodes` smallest
values are chosen.

"""
@kwdef struct GeneralizedRandomTesselated{I<:Integer} <: BONSampler
    num_nodes::I = _DEFAULT_NUM_NODES
end 


function _sample(
    sampler::GeneralizedRandomTesselated, 
    domain; 
    inclusion=nothing
)
    addresses = _get_addresses(sampler, domain)
    nodes = _pick_nodes(sampler, domain, addresses)
    return BiodiversityObservationNetwork(nodes, domain)
end 

"""
    _quadrant_fill!(mat)

Takes a matrix `mat` and splits it into quadrants randomly labeled one through four. 
"""
function _quadrant_fill!(mat)
    x, y = size(mat)
    a, b = x ÷ 2, y ÷ 2
    quad_ids = SB.shuffle([1, 2, 3, 4])
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
        mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]] .= _quadrant_fill!(mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]])
    end

    return mat
end

"""
    _get_address_length

Returns the number of digits in a [`GeneralizedRandomTessellatedStratified`](@ref)
address for a specific geometry, computed based on the raster dimensions. 
"""
_get_address_length(sampler, raster) = Int(ceil(max(log(2, size(raster,1)), log(2, size(raster,2)))))

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

"""
    _pick_nodes(sampler, raster, addresses)
"""
function _pick_nodes(sampler, raster, addresses)
    xbound, ybound = size(raster)
    sort_idx = sortperm([addresses[cidx] for cidx in eachindex(addresses)])
    cart_idxs = CartesianIndices(addresses)[sort_idx]
    
    num_selected = 0
    cursor = 1

    selected_nodes = []
    while num_selected < sampler.num_nodes && cursor <= length(sort_idx)
        candidate = cart_idxs[cursor]
        if candidate[1] <= xbound && candidate[2] <= ybound
            coord = CartesianIndex(candidate[1], candidate[2])
            
            if !ismasked(raster, coord)
                push!(selected_nodes, coord)
                num_selected += 1
            end
        end 
        cursor += 1
    end 
    return selected_nodes 
end
