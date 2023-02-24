module BONTestCubeSampling

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError CubeSampling(numpoints = 0)

# Must select fewer points than number of candidate points
@test_throws ArgumentError CubeSampling(numpoints = 20, pik = zeros(10), x = rand(0:4, 2, 10))

# Must have the same number of inclusion probabilites as auxillary variables
@test_throws DimensionMismatch CubeSampling(numpoints = 5, pik = zeros(15), x = rand(0:4, 2, 10))

# Correct subtype
@test typeof(CubeSampling(numpoints = 5, x = rand(0:4, 2, 10))) <: BONRefiner
@test typeof(CubeSampling(numpoints = 5, x = rand(0:4, 2, 10)))<: BONSampler

# Test with a random uncertainty matrix
N = 100
U = rand(20, 20)
pack = seed(BalancedAcceptance(; numpoints = N), U)
c = Vector{CartesianIndex}(undef, 15)
x_mat = rand(0:4, 2, N)
smpl = CubeSampling(numpoints = length(c), x = x_mat)

# Length and element type
@test length(first(refine(pack, smpl))) == smpl.numpoints
@test eltype(first(refine(first(pack), smpl, U))) == CartesianIndex

# Test with an existing coordinates vector
@test_throws DimensionMismatch refine!(
    c,
    first(pack),
    CubeSampling(; numpoints = length(c) - 1, x = x_mat),
    last(pack),
)

end 