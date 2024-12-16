module BONTestBalancedAcceptance

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws TooFewSites BalancedAcceptance(0, 1.0)

@test_throws TooManySites seed(BalancedAcceptance(; numsites = 26), rand(5, 5))

# Must have a positive alpha
@test_throws ArgumentError BalancedAcceptance(1, -0.01)

# Parametric constructor
@test typeof(BalancedAcceptance(2, 0.2f0)) == BalancedAcceptance{typeof(2), Float32}

# Correct subtype
@test typeof(BalancedAcceptance(2, 0.2)) <: BONSeeder
@test typeof(BalancedAcceptance(2, 0.2f0)) <: BONSampler

# Test with a random uncertainty matrix
U = rand(20, 20)
@test length(first(seed(BalancedAcceptance(; numsites = 10), U))) == 10
@test length(first(seed(BalancedAcceptance(; numsites = 20), U))) == 20
@test eltype(first(seed(BalancedAcceptance(; numsites = 10), U))) == CartesianIndex

# Test with an existing coordinates vector
c = Vector{CartesianIndex}(undef, 20)
@test_throws DimensionMismatch seed!(c, BalancedAcceptance(; numsites = length(c) - 1), U)

# Test the curried version
@test length(first(seed(BalancedAcceptance(; numsites = 12))(U))) == 12

# Test the curried allocating version
@test length(first(seed!(c, BalancedAcceptance(; numsites = length(c)))(U))) == length(c)
@test_throws DimensionMismatch seed!(c, BalancedAcceptance(; numsites = length(c) - 1))(U)

end
