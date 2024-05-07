module BONTestAdaptiveSpatial

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError AdaptiveSpatial(0)

# Correct subtype
@test typeof(AdaptiveSpatial(10)) <: BONRefiner
@test typeof(AdaptiveSpatial(10)) <: BONSampler

# Test with a random uncertainty matrix
U = rand(20, 20)
pack = seed(BalancedAcceptance(; numsites = 30), U)
c = Vector{CartesianIndex}(undef, 15)
smpl = AdaptiveSpatial(length(c))

# Length and element type
@test length(first(refine(pack, smpl))) == smpl.numsites
@test eltype(first(refine(first(pack), smpl, U))) == CartesianIndex

# Test with an existing coordinates vector
@test_throws DimensionMismatch refine!(
    c,
    first(pack),
    AdaptiveSpatial(; numsites = length(c) - 1),
    last(pack),
)

# Test the curried version
@test length(first(refine(AdaptiveSpatial(; numsites = 12))(pack...))) == 12

# Test the curried allocating version
@test length(first(refine!(c, AdaptiveSpatial(; numsites = length(c)))(pack...))) ==
      length(c)
@test_throws DimensionMismatch refine!(c, AdaptiveSpatial(; numsites = length(c) - 1))(
    pack...,
)

end