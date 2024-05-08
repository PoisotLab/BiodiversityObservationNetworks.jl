"""
    FractalTriad

A `BONSeeder` that generates `FractalTriad` designs
"""
Base.@kwdef struct FractalTriad{I <: Integer, F <: AbstractFloat} <: BONSeeder
    numsites::I = 81
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

maxsites(ft::FractalTriad) = maximum(ft.dims) * 10  # gets numerically unstable for large values because float coords belong to the the same cell in the raster, and arctan breaks  
function check_arguments(ft::FractalTriad)
    check(TooManySites, ft)
    ft.numsites > 9 || throw(TooFewSites("FractalTriad requires more than 9 or more sites"))
    ft.numsites % 9 == 0 ||
        throw(ArgumentError("FractalTriad requires number of sites to be a multiple of 9"))
    return
end

"""
    _outer_triangle

Takes a FractalTriad generator `ft` and returns the outermost triangle based on the size of the raster and the padding on the horizontal and vertical.
"""
function _outer_triangle(ft)
    x, y = ft.dims
    Δx, Δy =
        Int.([floor(0.5 * ft.horizontal_padding * x), floor(0.5 * ft.vertical_padding * y)])

    return CartesianIndex.([(Δx, Δy), (x ÷ 2, y - Δy), (x - Δx, Δy)])
end

"""
    _hexagon(A, B, C)

Takes a set of `CartesianIndex`'s that form the vertices of a triangle, and returns the `CartisianIndex`'s for each of the points that form the internal hexagon points

For input vertices A, B, C, `_hexagon` returns the points on the edges of the triangle below in the order `[d,e,f,g,h,i]`

                        B 

                  e           f


             d                      g

        A          i          h            C

--

After running `vcat(triangle, hex)`, the resulting indices form the 2-level triad with indices corresponding to points in the below manner:

                         2 


                  5           6


             4                      7

        1         9            8         3

  - 
γ = |AB|
χ = |BC|

θ = ⦤ BAC
α = ⦤ BCA

TODO: this always assumes |AC| is horizontal. This could be changed later.
"""
function _hexagon(A, B, C)
    γ = sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2) # left side length
    χ = sqrt((B[1] - C[1])^2 + (B[2] - C[2])^2) # right side length

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

"""
    _fill_triangle!(coords, traingle, count)

Takes a set of vertices of a triangle `triangle`, and fills the internal hexagon for those points.
"""
function _fill_triangle!(coords, triangle, count)
    start = count
    hex = _hexagon(triangle...)
    for i in eachindex(hex)
        coords[count] = CartesianIndex(hex[i])
        count += 1
        if count > length(coords)
            return coords[start:(count - 1)], count
        end
    end

    return coords[start:(count - 1)], count
end

"""
    _generate!(coords::Vector{CartesianIndex}, ft::FractalTriad)

Fills `coords` with a set of points generated using the `FractalTriad` generator `ft`.
"""
function _generate!(
    coords::Vector{CartesianIndex},
    ft::FractalTriad,
)
    base_triangle = _outer_triangle(ft)
    coords[1:3] .= base_triangle
    count = 4

    triangle = coords[1:3]
    hex, count = _fill_triangle!(coords, triangle, count)
    pack = vcat(triangle, hex)
    vert_idxs = [[5, 2, 6], [1, 4, 9], [8, 7, 3]]

    pack_stack = [pack]

    while length(pack_stack) > 0
        pack = popat!(pack_stack, 1)
        for idx in vert_idxs
            triangle = pack[idx]
            hex, count = _fill_triangle!(coords, triangle, count)
            if count > ft.numsites
                return coords
            end
            push!(pack_stack, vcat(triangle, hex))
        end
    end
    return coords
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "FractalTriad is correct subtype" begin
    @test FractalTriad <: BONSeeder
    @test FractalTriad <: BONSampler
end

@testitem "FractalTriad default constructor works" begin
    @test typeof(FractalTriad()) <: FractalTriad
end

@testitem "FractalTriad can change number of sites with keyword argument" begin
    ft = FractalTriad(; numsites = 18)
    @test typeof(ft) <: FractalTriad
    @test ft.numsites == 18

    ft = FractalTriad(; numsites = 27)
    @test typeof(ft) <: FractalTriad
    @test ft.numsites == 27

    @test_throws ArgumentError FractalTriad(numsites = 20)
end

@testitem "FractalTriad throws error when too few points as passed as an argument" begin
    @test_throws TooFewSites FractalTriad(; numsites = 9)
    @test_throws TooFewSites FractalTriad(; numsites = -1)
    @test_throws TooFewSites FractalTriad(; numsites = 0)
end
