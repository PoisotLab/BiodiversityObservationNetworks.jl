using Pkg
Pkg.activate(@__DIR__)

using BiodiversityObservationNetworks
using CairoMakie, GeoMakie
using NeutralLandscapes

import BiodiversityObservationNetworks as BONs
import BiodiversityObservationNetworks.SpeciesDistributionToolkit as SDT
import BiodiversityObservationNetworks.GeoInterface as GI
import BiodiversityObservationNetworks.GeometryOps as GO
import BiodiversityObservationNetworks.DelaunayTriangulation as DT
import BiodiversityObservationNetworks.SpeciesDistributionToolkit.SimpleSDMLayers.ArchGDAL as AGDAL


country_coda = "FRA"

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
ax, plt = bonplot(f[1,1], bon, _STATES, axistype=GeoAxis)
f


vor = voronoi(bon, _COUNTRY)

viridis = 
[Makie.ColorSchemes.viridis[i/length(vor)] for i in 1:length(vor)]


f = Figure()
ax = Axis(f[1,1])
for (i,v) in enumerate(vor)
    poly!(ax, v, color=(viridis[i], 0.7))
end
poly!(ax, _COUNTRY.geometry, color=(:blue, 0), strokewidth=1)
scatter!(ax, [n.coordinate for n in bon], color=:red)
f



f = Figure()
ax, plt = bonplot(f[1,1], bon, _STATES, axistype=GeoAxis)
poly!(p, strokewidth=2, color=(:white, 0))
f

tri = DT.triangulate([n.coordinate for n in bon])
vor = DT.voronoi(tri)

vor


f = Figure()
ax = GeoAxis(f[1,1])
triplot!(ax, res)
fs

bon = sample(AdaptiveHotspot(), bioclim[1])


