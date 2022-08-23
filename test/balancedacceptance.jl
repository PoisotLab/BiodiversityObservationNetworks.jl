module BONTestBalancedAcceptance

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError BalancedAcceptance(0, 1.0)

# Must have a positive alpha
@test_throws ArgumentError BalancedAcceptance(1, 0.0)

# Parametric constructor
@test typeof(BalancedAcceptance(2, 0.2f0)) == BalancedAcceptance{Int64,Float32}

end