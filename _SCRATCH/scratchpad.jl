using Pkg
Pkg.activate(@__DIR__)

using BiodiversityObservationNetworks
using CairoMakie, GeoMakie
using NeutralLandscapes

import BiodiversityObservationNetworks as BONs
import BiodiversityObservationNetworks.SpeciesDistributionToolkit as SDT
import BiodiversityObservationNetworks.GeoInterface as GI
import BiodiversityObservationNetworks.GeometryOps as GO


country_coda = "COL"

_COUNTRY = SDT.gadm(country_coda)
_STATES = SDT.gadm(country_coda, 1)


bioclim = SDT.SDMLayer[SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=i, SDT.boundingbox(_COUNTRY)...) for i in 1:19]
bioclim = RasterStack(SDT.SimpleSDMLayers.mask!(bioclim, _COUNTRY))


cornerplot(bioclim, )
``
bon = sample(SimpleRandom(50), bioclim)
bon = sample(SimpleRandom(5), _STATES)
bon = sample(SpatiallyStratified(100), _STATES)
bon = sample(Grid(), _COUNTRY)
bon = sample(SimpleRandom(100), _COUNTRY)
bon = sample(BalancedAcceptance(number_of_nodes=50), bioclim)
bon = sample(GeneralizedRandomTessellatedStratified(number_of_nodes=50), _COUNTRY)
bon = sample(GeneralizedRandomTessellatedStratified(number_of_nodes=50), _COUNTRY)
bon = sample(AdaptiveHotspot(), bioclim[1])

#bon = sample(UncertaintySampling(300), H)


# this should ensure simplerandom is applied on states-by-states basis by
# default 
bon = sample(MultistageSampler([BalancedAcceptance(20), SimpleRandom(10)]), _STATES)


bon = sample(CubeSampling(), bioclim)


bon = sample(MultistageSampler([BalancedAcceptance(250), CubeSampling()]), bioclim)

# NOTES: possible to imagine a situation where we want the number of BAS points
# in each state to be distributed by state area. There are several ways that
# this could be realized: 
# - SpatiallySpatified could take an argument that uses a specific sampler
#   within each spatial strata, rather than always SRS.
bon = sample(BalancedAcceptance(200), bioclim)


f = Figure(size=(500, 500))
bonplot(f[1,1], bon, _STATES, axistype=GeoAxis)
f





# How to measure the distance between a BON and the whole env space in layer
# stack?
# Marginals instead of MvNormal bc everything gets fucky 
# But can penalize w/ weighted distance between Covariance matrices

function _js_thing(rs::RasterStack, bon::BiodiversityObservationNetwork)
    function _js(P,Q) 
        M = BONs.Distributions.MixtureModel([P,Q], [0.5,0.5])
        div = 0.5*BONs.Distributions.kldivergence(P,M) + 0.5*BONs.Distributions.kldivergence(Q,M)
        return sqrt(div / log(2))    
    end
    _, Xfull = BONs.features(rs)
    Xsampled = rs[bon]

    nlayers = length(rs)
    
    Σ₁, Σ₂ = zeros(nlayers, nlayers), zeros(nlayers, nlayers)

    𝓛_js = 0.
    for i in axes(Xfull,1)
        𝓝₁ = BONs.Distributions.fit(BONs.Distributions.Normal, Xfull[i,:])
        𝓝₂ = BONs.Distributions.fit(BONs.Distributions.Normal, Xsampled[i,:])
        𝓛_js += _js(𝓝₁, 𝓝₂) 
        for j in i+1:nlayers
            Σ₁[i,j] = BONs.Distributions.cov(Xfull[:,i],Xfull[:,j])
            Σ₁[j,i] = Σ₁[i,j]

            Σ₂[i,j] = BONs.Distributions.cov(Xsampled[:,i],Xsampled[:,j])
            Σ₂[j,i] = Σ₂[i,j]
        end 
    end
    𝓛_covariance = sqrt(sum((Σ₁ .- Σ₂).^2))

    𝓛_js, 𝓛_covariance
end 


n_reps = 500
n_nodes = [2^i for i in 4:9]
begin 
    f = Figure()
    ax = Axis(f[1,1])
    for n in n_nodes
        density!(ax, [_js_thing(bioclim, sample(SimpleRandom(n), bioclim))[1] for _ in 1:n_reps])
    end 
    f
end 