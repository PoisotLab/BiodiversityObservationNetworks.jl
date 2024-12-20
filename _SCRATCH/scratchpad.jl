using Pkg
Pkg.activate(@__DIR__)

using BiodiversityObservationNetworks
using CairoMakie, GeoMakie


import BiodiversityObservationNetworks as BONs
import BiodiversityObservationNetworks.SpeciesDistributionToolkit as SDT
import BiodiversityObservationNetworks.GeoInterface as GI
import BiodiversityObservationNetworks.GeometryOps as GO


col = SDT.gadm("COL")
col_states = SDT.gadm("COL", 1)

bioclim = RasterStack([SDT.mask!(SDT.SDMLayer(SDT.RasterData(SDT.WorldClim2, SDT.BioClim); layer=i, SDT.boundingbox(col)...), col) for i in 1:19])

bon = sample(SimpleRandom(50), bioclim)
bon = sample(SimpleRandom(10), col_states)
bon = sample(SpatiallyStratified(100), col_states)
bon = BONs.sample(Grid(), col)
bon = BONs.sample(KMeans(75), bioclim)


cornerplot(bioclim)



begin 
    f = Figure(size=(900,900))
    ga = GeoAxis(f[1,1])
    for (i,st) in enumerate(convert(Polygon, col_states))
    poly!(ga, st.geometry, strokewidth=2, color=Symbol("grey", string(20+2*i)))
    scatter!(ga, [node[1] for node in bon], [node[2] for node in bon], color=(:red))
    end
    f
end 

