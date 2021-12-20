
mutable struct BON
    numobservatories
    coordinates 
end


abstract type SpatialSampler end 

struct FractalTriad <: SpatialSampler end
struct SpatialSimulatedAnnealing <: SpatialSampler end
struct BalancedAcceptanceSampling <: SpatialSampler end