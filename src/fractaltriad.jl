Base.@kwdef struct FractalTriad{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numsites::I = 30
    horizontal_padding::F = 0.1
    vertical_padding::F = 0.1
    dims::Tuple{I, I} = (100, 100)
    function FractalTriad(numsites, horizontal_padding, vertical_padding, dims)
        ft = new{typeof(numsites), typeof(horizontal_padding)}(
            numsites,
            horizontal_padding,
            vertical_padding,
            dims,
        )
        check_arguments(ft)
        return ft
    end
end

maxsites(ft::FractalTriad) = prod(ft.dims)

function check_arguments(ft::FractalTriad)
    check(TooManySites, ft)
    return check(TooFewSites, ft)
end

function _outer_triangle(ft)
    x, y = ft.dims
    Δx, Δy =
        Int.([floor(0.5 * ft.horizontal_padding * x), floor(0.5 * ft.vertical_padding * y)])

    return CartesianIndex.([(Δx, Δy), (x ÷ 2, y - Δy), (x - Δx, Δy)])
end

"""
This is the layout

                        B 

                  e           f


             d                      g

        A          i          h            C

---

                        2 

                  5           6


             4                      7

        1         9            8         3


θ = ∠ BAC

d = γ/3 cos(θ), γ/3 sin(θ)
e = γ 

"""

function _hexagon(outside)
    A, B, C = outside

    γ = sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2)
    χ = sqrt((B[1] - C[1])^2 + (B[2] - C[2])^2)

    θ = atan((B[2] - A[2]) / (B[1] - A[1]))
    α = atan((B[2] - C[2]) / (C[1] - B[1]))

    d = (A[1] + (γ * cos(θ) / 3), A[2] + γ * sin(θ) / 3)
    e = (A[1] + 2γ * cos(θ) / 3, A[2] + 2γ * sin(θ) / 3)
    f = (A[1] + (C[1] - A[1]) - (2χ * cos(α) / 3), A[2] + 2χ * sin(α) / 3)
    g = (A[1] + (C[1] - A[1] - (χ * cos(α) / 3)), A[2] + χ * sin(α) / 3)
    h = (A[1] + 2(C[1] - A[1]) / 3, A[2])
    i = (A[1] + (C[1] - A[1]) / 3, A[2])
    return [CartesianIndex(Int.([floor(x[1]), floor(x[2])])...) for x in [d, e, f, g, h, i]]
end

function _fill_triangle(coords, triangle, count)
    start = count
    hex = _hexagon(triangle)
    for i in eachindex(hex)
        coords[count] = CartesianIndex(hex[i])
        count += 1
        if count > length(coords)
            return coords[start:(count - 1)], count
        end
    end

    return coords[start:(count - 1)], count
end

function _generate!(
    coords::Vector{CartesianIndex},
    ft::FractalTriad,
)
    base_triangle = _outer_triangle(ft)
    coords[1:3] .= base_triangle

    count = 4

    triangle = coords[1:3]
    hex, count = _fill_triangle(coords, triangle, count)

    pack = vcat(triangle, hex)
    vert_idxs = [[5, 2, 6], [1, 4, 9], [8, 7, 3]]

    pack_stack = [pack]

    while length(pack_stack) > 0
        pack = popat!(pack_stack, 1)
        for idx in vert_idxs
            triangle = pack[idx]
            hex, count = _fill_triangle(coords, triangle, count)
            if count > ft.numsites
                return coords
            end
            push!(pack_stack, vcat(triangle, hex))
        end
    end
    return coords
end
