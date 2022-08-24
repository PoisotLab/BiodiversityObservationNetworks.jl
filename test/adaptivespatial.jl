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
pool = seed(BalancedAcceptance(numpoints=30), U)
c = Vector{CartesianIndex}(undef, 15)
smpl = AdaptiveSpatial(length(c))

# Length and element type
@test length(refine(pool, smpl, U)) == smpl.numpoints
@test eltype(refine(pool, smpl, U)) == CartesianIndex

# Test with an existing coordinates vector
@test_throws DimensionMismatch refine!(c, pool, AdaptiveSpatial(numpoints=length(c) - 1), U)

# Test the curried version
@test length(refine(AdaptiveSpatial(numpoints=12))(pool, U)) == 12

# Test the curried allocating version
@test length(refine!(c, AdaptiveSpatial(numpoints=length(c)))(pool, U)) == length(c)
@test_throws DimensionMismatch refine!(c, AdaptiveSpatial(numpoints=length(c) - 1))(pool, U)

end