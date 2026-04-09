"""
    GRTS <: BONSampler

Generalized Random Tessellation Stratified (GRTS) sampling.

GRTS produces spatially balanced samples by recursively partitioning the
domain into quadrants, assigning each cell a hierarchical address. The recursive
tessellation ensures spatial spread without requiring distance computations.

Originally proposed by Stevens & Olsen (2004) (CITE TODO).

# Fields
- `n::Int`: number of sites to select (default 50)
"""
@kwdef struct GRTS <: BONSampler
    n::Int = 50
end

guarantees_exact_n(::GRTS) = true

"""
    _grts_grid_size

Determine grid dimensions for the GRTS tessellation. Uses metadata if
available, otherwise infers from coordinate extent.
"""
function _grts_grid_size(cpool::CandidatePool)
    xmax, ymax = maximum([x[1] for x in cpool.keys]), maximum([x[2] for x in cpool.keys])
    return (xmax, ymax)
end

"""
    _get_address_length

Returns the number of digits in a [`GeneralizedRandomTessellatedStratified`](@ref)
address for a specific geometry, computed based on the raster dimensions. 
"""
_get_address_length(sz) = ceil(Int, max(log2(sz[1]), log2(sz[2])))

"""
    _quadrant_split!(mat, grid_size)

Splits a matrix `mat` into nested quadrants, where the side-length of a
submatrix to be split into quadrants is given by `grid_size`.  
"""
function _quadrant_split!(rng, mat, grid_size)
    x, y = size(mat)
    num_x_grids, num_y_grids = x ÷ grid_size, y ÷ grid_size

    for i in 0:(num_x_grids - 1), j in 0:(num_y_grids - 1)
        a, b = i * grid_size, j * grid_size
        bounds = (a + 1, b + 1), (a + grid_size, b + grid_size)
        mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]] .= _quadrant_fill!(rng, mat[bounds[1][1]:bounds[2][1], bounds[1][2]:bounds[2][2]])
    end

    return mat
end

"""
    _quadrant_fill!(mat)

Takes a matrix `mat` and splits it into quadrants randomly labeled one through four. 
"""
function _quadrant_fill!(rng, mat)
    x, y = size(mat)
    a, b = x ÷ 2, y ÷ 2
    quad_ids = StatsBase.shuffle(rng, [1, 2, 3, 4])
    mat[begin:a, begin:b] .= quad_ids[1]
    mat[(a + 1):end, begin:b] .= quad_ids[2]
    mat[begin:a, (b + 1):end] .= quad_ids[3]
    mat[(a + 1):end, (b + 1):end] .= quad_ids[4]
    return mat
end

"""
    _get_addresses(sz)

Gets an address matrix of size `sz`. 
"""
function _get_addresses(rng, sz)
    address_length = _get_address_length(sz)
    grid_sizes = reverse([2^i for i in 1:address_length])
    grids = [zeros(Int, 2^address_length, 2^address_length) for _ in grid_sizes]

    map(
        i -> _quadrant_split!(rng, grids[i], grid_sizes[i]),
        eachindex(grids),
    )
    addresses = sum([10^(i - 1) .* ag for (i, ag) in enumerate(grids)])
    return addresses
end



"""
    _pick_nodes(sampler, raster, addresses)
"""
function _pick_nodes(sampler, cpool, addresses, grid_size)
    mask, _, key_index =_construct_explicit_mask_and_inclusion(cpool, grid_size)

    sort_idx = sortperm([addresses[cidx] for cidx in eachindex(addresses)])
    cart_idxs = CartesianIndices(addresses)[sort_idx]
    
    num_selected = 0
    cursor = 1

    selected = []
    while num_selected < sampler.n && cursor <= length(sort_idx)
        candidate = cart_idxs[cursor]
        if candidate[1] <= grid_size[1] && candidate[2] <= grid_size[2]
            coord = CartesianIndex(candidate[1], candidate[2])
            if mask[coord] 
                push!(selected, key_index[coord])
                num_selected += 1
            end
        end 
        cursor += 1
    end 
    return selected 
end

function _sample(rng::AbstractRNG, sampler::GRTS, cpool::CandidatePool)
    grid_size = _grts_grid_size(cpool)
    addresses = _get_addresses(rng, grid_size)
    _pick_nodes(sampler, cpool, addresses, grid_size)
end 

@testitem "GRTS works" begin
    bon = sample(GRTS(), rand(30,20))
    @test bon isa BiodiversityObservationNetwork
end
