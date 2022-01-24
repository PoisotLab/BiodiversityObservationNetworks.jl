
@kwdef mutable struct BiodiversityObservationNetwork{IT<:Integer,M<:AbstractMatrix}
    numobservatories::IT = 50
    coordinates::M = missing
end


abstract type SpatialSampler end 

struct SpatialSimulatedAnnealing <: SpatialSampler end