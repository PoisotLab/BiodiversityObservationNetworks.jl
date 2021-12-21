
mutable struct BON
    numobservatories
    coordinates 
end


abstract type SpatialSampler end 

struct SpatialSimulatedAnnealing <: SpatialSampler end
struct BalancedAcceptanceSampling <: SpatialSampler end