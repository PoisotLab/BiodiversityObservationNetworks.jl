module BONTestUniqueness

using BiodiversityObservationNetworks
using Test

# Must have one point or more
@test_throws ArgumentError Uniqueness(0, zeros(5,5,5))
@test_throws ArgumentError Uniqueness(0, zeros(5,5))

# Parametric constructor
@test typeof(Uniqueness(Int32(2), zeros(Float32, 5,5,5))) == Uniqueness{Int32,Float32}

# Correct subtype
@test typeof(Uniqueness(2, zeros(5,5,5))) <: BONRefiner
@test typeof(Uniqueness(2, zeros(5,5,5))) <: BONSampler

# Test with a random uncertainty matrix
layers = rand(20, 20, 5)
np = 25
pack = seed(BalancedAcceptance(; numpoints = np), rand(20,20))

@test length(first(refine(pack, Uniqueness(; numpoints = 10, layers=layers)))) == 10
@test length(first(refine(pack, Uniqueness(; numpoints = 20, layers=layers)))) == 20
@test eltype(first(refine(pack, Uniqueness(; numpoints = 10, layers=layers)))) == CartesianIndex

# Test with an existing coordinates vector
c = Vector{CartesianIndex}(undef, np)
@test_throws DimensionMismatch refine!(
    c,
    first(pack),
    Uniqueness(; numpoints = length(c) - 1, layers=layers),
    last(pack),
)
# Test the curried allocating version
@test length(first(refine(Uniqueness(; numpoints = 12, layers=layers))(pack...))) == 12

# Test the curried allocating version
@test length(first(refine!(c, Uniqueness(; numpoints = length(c), layers=layers))(pack...))) ==
      length(c)
@test_throws DimensionMismatch refine!(c, Uniqueness(; numpoints = length(c) - 1, layers=layers))(
    pack...,
)
end