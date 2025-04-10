
# Include unit j at step j with probability πⱼ. 
# If included, set Iᵢ = 1, else set Iᵢ = 0 

# Update rule for unit i that has not yet been considered after considering unit j
# πᵢ⁽ʲ⁾ = πᵢ⁽ʲ⁻¹⁾ - ( Iⱼ - πⱼ⁽ʲ⁻¹⁾ ) wⱼ⁽ⁱ⁾

# The weight given to the unit currently being considered, j, by each unit not
# yet considered, i, denoted wⱼ⁽ⁱ⁾ is bounded by (via the update rule)
#  -min( (1-πᵢ⁽ʲ⁻¹⁾) / (1-πⱼ⁽ʲ⁻¹⁾)  , πᵢ⁽ʲ⁻¹⁾ / πⱼ⁽ʲ⁻¹⁾  ) <= wⱼ⁽ⁱ⁾ 
# <= min( πᵢ⁽ʲ⁻¹⁾ / (1 - πⱼ⁽ʲ⁻¹⁾)  ,  (1 - πᵢ⁽ʲ⁻¹⁾) / πⱼ⁽ʲ⁻¹⁾ ) 

# Maximum weighting strategy: after deciding whether or not to include unit j, choose maximum possible
# weight for each subsequent unit i > j, but in order of how close they are to j

""" 
    SpatiallyCorrelatedPoisson

Spatially Correlated Poisson Sampling [Grafstrom2012SpaCor](@cite)

"""
Base.@kwdef struct SpatiallyCorrelatedPoisson{I<:Integer} <: BONSampler
    number_of_nodes::I = _DEFAULT_NUM_NODES
end 

_valid_geometries(::SpatiallyCorrelatedPoisson) = (BiodiversityObservationNetwork)

_max_weight(Π, i, j) = min(Π[i]/(1-Π[j]), (1-Π[i])/(Π[j]))

function _sample(
    sampler::SpatiallyCorrelatedPoisson, 
    bon::BiodiversityObservationNetwork,
    Π = fill(sampler.number_of_nodes/size(bon), size(bon))
) 
    inclusion_indicator = zeros(Bool, size(bon))
    Πⱼ₋₁ = deepcopy(Π)
    weights = zeros(size(bon))

    dist_mat = _get_distance_matrix(bon)

    closest_idx_mat = vcat([sortperm(r)[2:end] for r in eachrow(dist_mat)]'...)

    for j in 1:size(bon)
        inclusion_indicator[j] = rand() < Π[j]
        Iⱼ = inclusion_indicator[j]

        weights .= 0
        Πⱼ₋₁ .= Π
            
        for i in closest_idx_mat[j,:]
            if i > j 
                weights[i] = _max_weight(Πⱼ₋₁, i, j)
            end
        end 
        
        weights ./= sum(weights)

        for i in closest_idx_mat[j,:]
            if i > j 
                Π[i] = Πⱼ₋₁[i] - (Iⱼ - Πⱼ₋₁[j]) * weights[i]
            end 
        end 
    end    
    


    return BiodiversityObservationNetwork(bon[inclusion_indicator])
end 

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use SCPS with default constructor on a BON" begin
    #polygon = openstreetmap("Colombia")
    layer = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(rand(150, 150))
    candidate_bon = sample(SimpleRandom(300), layer)

    scps = SpatiallyCorrelatedPoisson()
    bon = sample(scps, candidate_bon)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == scps.number_of_nodes
end