module BONTestAdaptiveSpatialSampling

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError AdaptiveSpatialSampling(0)

# Correct subtype
@test typeof(AdaptiveSpatialSampling(10)) <: BONRefiner
@test typeof(AdaptiveSpatialSampling(10)) <: BONSampler

# Test with a random uncertainty matrix
U = rand(20, 20)
pool = seed(BalancedAcceptance(numpoints=30), U)
c = Vector{CartesianIndex}(undef, 15)
smpl = AdaptiveSpatialSampling(length(c))

# Length and element type
@test length(refine(pool, smpl, U)) == smpl.numpoints
@test eltype(refine(pool, smpl, U)) == CartesianIndex

# Test with an existing coordinates vector
@test_throws DimensionMismatch refine!(c, pool, AdaptiveSpatialSampling(numpoints=length(c) - 1), U)

# Test the curried version
@test length(refine(AdaptiveSpatialSampling(numpoints=12))(pool, U)) == 12

# Test the curried allocating version
@test length(refine!(c, BalancedAcceptance(numpoints=length(c)))(U)) == length(c)
@test_throws DimensionMismatch refine!(c, BalancedAcceptance(numpoints=length(c) - 1))(U)

end