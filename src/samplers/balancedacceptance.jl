"""
    BalancedAcceptance

`BalancedAcceptance` is a type of [`BONSampler`](@ref) for generating
[`BiodiversityObservationNetwork`](@ref)s with spatial spreading.

*Arguments*:
- `number_of_nodes`: the number of sites to select
- `grid_size`: if being used on a polygon, the dimensions of the grid used to
  cover the extent. Balanced Acceptance sampling uses discrete Cartesian indices
"""
Base.@kwdef struct BalancedAcceptance{I<:Integer} <: BONSampler
    number_of_nodes::I = 100
    grid_size::Tuple{I,I} = (250, 250)
end 

_valid_geometries(::BalancedAcceptance) = (Polygon, Raster, Vector{Polygon}, RasterStack)

function _sample(sampler::BalancedAcceptance, raster::T) where T<:Union{Polygon,Raster,Vector{<:Polygon}}
    _balanced_acceptance(sampler, raster)
end 

function _sample(sampler::BalancedAcceptance, layers::RasterStack)
    _balanced_acceptance(sampler, first(layers))
end 

_get_easting_and_northing(::BalancedAcceptance, raster::Raster) = SDT.eastings(raster), SDT.northings(raster)
_get_easting_and_northing(sampler::BalancedAcceptance, layers::RasterStack) = _get_easting_and_northing(sampler, first(layers))

function _get_easting_and_northing(sampler::BalancedAcceptance, polygon::Polygon)
    x, y = GI.extent(polygon)
    grid_size = sampler.grid_size
    Δx = (x[2]-x[1])/grid_size[1]
    Δy = (y[2]-y[1])/grid_size[2]

    Es = [x[1] + i*Δx for i in 1:grid_size[1]]
    Ns = [y[1] + i*Δy for i in 1:grid_size[2]]

    return Es, Ns
end 

_check_candidate(Es, Ns, candidate, polygon::Polygon) = GeometryOps.contains(polygon, (Es[candidate[1]], Ns[candidate[2]]))

function _check_candidate(Es, Ns, candidate, raster::Raster)
    val = raster.raster[candidate[2], candidate[1]]
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



#=
"""
    BalancedAcceptance

A `BONSeeder` that uses Balanced-Acceptance Sampling (Van-dem-Bates et al. 2017
https://doi.org/10.1111/2041-210X.13003)
"""
Base.@kwdef struct BalancedAcceptance{I<:Integer} <: BONSampler
    numsites::I = 30
    dims::Tuple{I,I} = (50,50)
    function BalancedAcceptance(numsites::I, dims::Tuple{J,J}) where {I<:Integer, J<:Integer}
        bas = new{I}(numsites, dims)
        check_arguments(bas) 
        return bas
    end
end

_default_pool(bas::BalancedAcceptance) = pool(bas.dims)
BalancedAcceptance(M::Matrix{T}; numsites = 30) where T = BalancedAcceptance(numsites, size(M))
BalancedAcceptance(l::Layer; numsites = 30) = BalancedAcceptance(numsites, size(l))

maxsites(bas::BalancedAcceptance) = prod(bas.dims)

function check_arguments(bas::BalancedAcceptance)
    check(TooFewSites, bas)
    check(TooManySites, bas)
    return nothing
end

function _sample!(
    selected::S,
    candidates::C,
    ba::BalancedAcceptance
) where {S<:Sites,C<:Sites}
    seed = rand(Int32.(1e0:1e7), 2)
    n = numsites(ba)
    x,y = ba.dims

    candidate_mask = zeros(Bool, x, y)
    candidate_mask[candidates.coordinates] .= 1

    # This is sequentially adding points, needs to check if that value is masked
    # at each step and skip if so  
    exp_needed = 10 * Int(ceil((length(candidates)/(x*y)) .* n))

    ct = 1
    for ptct in 1:exp_needed
        i, j = haltonvalue(seed[1] + ptct, 2), haltonvalue(seed[2] + ptct, 3)
        proposal = CartesianIndex(convert.(Int32, [ceil(x * i), ceil(y * j)])...)
        if ct > n 
            break
        end 
        if candidate_mask[proposal]
            selected[ct] = proposal
            ct += 1
        end 
    end
    return selected
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "BalancedAcceptance default constructor works" begin
    @test typeof(BalancedAcceptance()) <: BalancedAcceptance
end

@testitem "BalancedAcceptance requires positive number of sites" begin
    @test_throws TooFewSites BalancedAcceptance(numsites = 1)
    @test_throws TooFewSites BalancedAcceptance(numsites = 0)
    @test_throws TooFewSites BalancedAcceptance(numsites = -1)
end

@testitem "BalancedAcceptance can't be run with too many sites" begin
    numpts, numcandidates = 26, 25
    @test numpts > numcandidates   # who watches the watchmen?
    dims = Int32.(floor.((sqrt(numcandidates), sqrt(numcandidates))))
    @test_throws TooManySites BalancedAcceptance(numpts, dims)
end

@testitem "BalancedAcceptance can generate points" begin
    bas = BalancedAcceptance()
    coords = sample(bas)

    @test typeof(coords) <: Sites
end


@testitem "BalancedAcceptance can take number of points as keyword argument" begin
    N = 40
    bas = BalancedAcceptance(; numsites = N)
    @test bas.numsites == N
end

=#