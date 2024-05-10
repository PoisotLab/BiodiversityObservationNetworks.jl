"""
    GeneralizedRandomTessellatedStratified

@Olsen
"""
@kwdef struct GeneralizedRandomTessellatedStratified <: BONSeeder
    numsites = 50
    dims = (100, 100)
end

maxsites(grts::GeneralizedRandomTessellatedStratified) = prod(dims)

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

    return coords .= CartesianIndices(code_numbers)[sort_idx][1:(grts.numsites)]
end
