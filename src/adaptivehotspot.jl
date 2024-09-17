"""
    AdaptiveHotspot

- **numsites**: an integer specifying the number of points to use (default is 30)
- **pool**: the sites that could potentially be picked 
- **uncertainty**: a `Layer` specifying the current uncertainty at each site
"""
Base.@kwdef mutable struct AdaptiveHotspot{T <: Integer, F <: AbstractFloat} <: BONSampler
    numsites::T = 30
    uncertainty::Layer = Layer(rand(50, 50))
    function AdaptiveHotspot(numsites, uncertainty)
        as = new{typeof(numsites), typeof(uncertainty[begin])}(numsites, uncertainty)
        check_arguments(as)
        return as
    end
end

AdaptiveHotspot(uncertainty::Matrix{T}; numsites = 30) where T = AdaptiveHotspot(numsites, uncertainty)

maxsites(as::AdaptiveHotspot) = prod(size(as.uncertainty))
function check_arguments(as::AdaptiveHotspot)
    check(TooFewSites, as)
    check(TooManySites, as)
    return nothing
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::AdaptiveHotspot,
)
    uncertainty = sampler.uncertainty
    # Distance matrix (inlined)
    d = zeros(Float64, Int((sampler.numsites * (sampler.numsites - 1)) / 2))

    # Start with the point with maximum entropy
    imax = last(findmax([uncertainty[i] for i in pool]))
    # Add it to the stack
    coords[1] = popat!(pool, imax)

    best_score = 0.0
    best_s = 1

    for i in 2:(sampler.numsites)
        for (ci, cs) in enumerate(pool)
            coords[i] = cs
            # Distance update
            start_from = Int((i - 1) * (i - 2) / 2) + 1
            end_at = start_from + Int(i - 2)
            d_positions = start_from:end_at
            for ti in 1:(i - 1)
                d[d_positions[ti]] = _D(cs, coords[ti])
            end
            # Get the score
            score = uncertainty[cs] + sqrt(log(i)) * _h(d[1:end_at], 1.0, 0.5)
            if score > best_score
                best_score = score
                best_s = ci
            end
        end
        coords[i] = popat!(pool, best_s)
    end
    return coords
end

function _matérn(d, ρ, ν)
    # This is the version from the supp mat
    # ν = 0.5 to have the exponential version
    return 1.0 * (2.0^(1.0 - ν)) / gamma(ν) *
           (sqrt(2ν) * d / ρ)^ν *
           besselk(ν, sqrt(2ν) * d / ρ)
end

function _h(d, ρ, ν)
    K = [_matérn(i, ρ, ν) for i in d]
    return (0.5 * log(2 * π * ℯ)^length(d)) * sum(K)
end

function _D(a1::T, a2::T) where {T <: CartesianIndex{2}}
    x1, y1 = a1.I
    x2, y2 = a2.I
    return sqrt((x1 - x2)^2.0 + (y1 - y2)^2.0)
end

# ====================================================
#
#   Tests
#
# =====================================================

@testitem "AdaptiveHotspot default constructor works" begin
    @test typeof(AdaptiveHotspot()) <: AdaptiveHotspot
end

@testitem "AdaptiveHotspot has right subtypes" begin
    @test AdaptiveHotspot <: BONSampler
end

@testitem "AdaptiveHotspot requires positive number of sites" begin
    @test_throws TooFewSites AdaptiveHotspot(numsites = 1)
    @test_throws TooFewSites AdaptiveHotspot(numsites = 0)
    @test_throws TooFewSites AdaptiveHotspot(numsites = -1)
end
