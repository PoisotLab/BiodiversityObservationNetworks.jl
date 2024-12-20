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