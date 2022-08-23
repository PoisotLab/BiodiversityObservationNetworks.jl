module BONTestBalancedAcceptance

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError BalancedAcceptance(0, 1.0)

# Must have a positive alpha
@test_throws ArgumentError BalancedAcceptance(1, 0.0)

# Parametric constructor
@test typeof(BalancedAcceptance(2, 0.2f0)) == BalancedAcceptance{Int64,Float32}

# Correct subtype
@test typeof(BalancedAcceptance(2, 0.2f0)) <: BONSeeder
@test typeof(BalancedAcceptance(2, 0.2f0)) <: BONSampler

# Test with a random uncertainty matrix
U = rand(20, 20)
@test length(seed(BalancedAcceptance(numpoints=10), U)) == 10
@test length(seed(BalancedAcceptance(numpoints=20), U)) == 20
@test eltype(seed(BalancedAcceptance(numpoints=10), U)) == CartesianIndex

# Test with an existing coordinates vector
c = Vector{CartesianIndex}(undef, 20)
@test_throws DimensionMismatch seed!(c, BalancedAcceptance(numpoints=length(c) - 1), U)

# Test the curried version
@test length(seed(BalancedAcceptance(numpoints=12))(U)) == 12

# Test the curried allocating version
@test length(seed!(c, BalancedAcceptance(numpoints=length(c)))(U)) == length(c)
@test_throws DimensionMismatch seed!(c, BalancedAcceptance(numpoints=length(c)-1))(U)

end