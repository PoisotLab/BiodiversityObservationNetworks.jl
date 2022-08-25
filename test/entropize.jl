module BONTestEntropize

using BiodiversityObservationNetworks
using Test

# Test with a random uncertainty matrix
U = rand(20, 20)
E = entropize(U)

@test minimum(E) == zero(eltype(E))
@test maximum(E) == one(eltype(E))

entropize!(E, U)

@test minimum(E) == zero(eltype(E))
@test maximum(E) == one(eltype(E))

end