
# Measuring Spatial Balance {#Measuring-Spatial-Balance}

BiodiversityObservationMethods contains a few utilities for quantifying the spatial balance of a set of nodes in a [`BiodiversityObservationNetwork`](/reference/api#BiodiversityObservationNetworks.BiodiversityObservationNetwork).

These methods are all used by calling the [`spatialbalance`](/reference/api#BiodiversityObservationNetworks.spatialbalance-Tuple{MoransI,%20BiodiversityObservationNetwork,%20RasterDomain}) method
<details class='jldocstring custom-block' open>
<summary><a id='BiodiversityObservationNetworks.spatialbalance-reference-utilities-spatialbalance' href='#BiodiversityObservationNetworks.spatialbalance-reference-utilities-spatialbalance'><span class="jlbinding">BiodiversityObservationNetworks.spatialbalance</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
spatialbalance
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/utilities/spatialbalance.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>



```julia
spatialbalance(::MoransI, raster::SDMLayer, bon::BiodiversityObservationNetwork)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/utilities/spatialbalance.jl#L58-L60" target="_blank" rel="noreferrer">source</a></Badge>



```julia
spatialbalance(::Type{VoronoiVariance}, bon::BiodiversityObservationNetwork, geom)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/blob/cf81ee6a79734664b33fd1821e1d682f012a7cc2/src/utilities/spatialbalance.jl#L133-L135" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Voronoi Variance {#Voronoi-Variance}

## Adjusted Moran&#39;s I {#Adjusted-Moran's-I}
