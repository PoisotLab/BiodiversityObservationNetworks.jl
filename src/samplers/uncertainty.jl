"""
    UncertaintySampling

Active-learning method. 
"""
Base.@kwdef struct UncertaintySampling{I<:Integer} <: BONSampler
    number_of_nodes::I = _DEFAULT_NUM_NODES
end 
_valid_geometries(::UncertaintySampling) = (SDMLayer)

function _is_uncertainty_layer(raster::SDMLayer) 
    Is = findall(raster.indices)
    sum(raster[Is] .> 0) > 0 && length(findall(x->x<0, raster[Is])) == 0
end

"""
    _sample(sampler::UncertaintySampling, geometry)

Internal dispatch for sampling using
[`UncertaintySampling`](@ref) on an `SDMLayer` containing uncertainty values.
"""
function _sample(sampler::UncertaintySampling, uncertainty::SDMLayer) 
    _is_uncertainty_layer(uncertainty) || throw(ArgumentError("Supplied raster is not a valid uncertainty layer."))

    N = sampler.number_of_nodes

    valid_idxs = findall(uncertainty.indices)
    
    weights = StatsBase.Weights(uncertainty[valid_idxs])
        

    cart_idxs = StatsBase.sample(valid_idxs, weights, N)

    Es, Ns = SDT.eastings(uncertainty), SDT.northings(uncertainty)
    BiodiversityObservationNetwork([Node(Es[c[2]], Ns[c[1]]) for c in cart_idxs])
end 

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use UncertaintySampling with default constructor on a Raster" begin
    uncert = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(rand(50,30))
    us = UncertaintySampling()
    bon = sample(us, uncert)
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == us.number_of_nodes
end

