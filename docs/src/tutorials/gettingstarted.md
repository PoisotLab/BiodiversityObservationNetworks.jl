# Getting Started with BiodiversityObservationNetworks.jl

In this tutorial, we will cover the basics of how to use BiodiversityObservationNetworks.jl. We'll start by loading the package.

```@ansi 1
using BiodiversityObservationNetworks
```

## "Hello World" in BiodiversityObservationNetworks.jl

The primary use of BiodiversityObservationNetworks is for generating [`BiodiversityObservationNetwork`](@ref) using a variety of point-selection algorithms. Generating networks, regardless of the specific algorithm chosen, is done using the [`sample`](@ref) method. The simplest point-selection algorithm is [`SimpleRandom`](@ref), where each location in space (what is meant be location? More on this in the next section) has an equal probability of inclusion.

```@example 1
bon = sample(SimpleRandom())
```

By default, [`SimpleRandom`](@ref) (and every other sampler), chooses 50 points. Without any other arguments, [`sample`](@ref) chooses points from a raster than covers the entire globe, with each pixel representing a 1˚ by 1˚ region.

When a version of the [Makie](https://docs.makie.org/v0.22/) package for data visualization is loaded, we can use built-in functions to visualize the network. Let's use [`GeoMakie`](https://geo.makie.org/v0.7.9/ ), which specifically enables plotting of geographic plots within Makie, with the `CairoMakie` backend. (Learn mroe about Makie backends [here](https://docs.makie.org/stable/explanations/backends/backends#What-is-a-backend))


```@example 1
using CairoMakie, GeoMakie
```

Now that these packages are loaded, we can use the [`bonplot`](@ref) method to visualize the generated [`BiodiversityObservationNetwork`](@ref). We pass `GeoAxis` to the keyword argument `axistype` to ensure [`bonplot`](@ref) uses `GeoMakie`.


```@example 1
bonplot(bon, axistype = GeoAxis)
```

We can adjust the number of points to generate by passing an integer directly to [`SimpleRandom`](@ref), i.e.

```@example 1
srs = SimpleRandom(150)
```

Alternatively, we can use the the `number_of_nodes` keyword argument

```@example 1
srs = SimpleRandom(number_of_nodes=150)
```

Both of these methods for adjusting the number of nodes is supported for all sampling algorithms.

Let's sample and visualize a [`BiodiversityObservationNetwork`](@ref) with more points

```
bon = sample(srs)
bonplot(bon, axistype=GeoAxis)
```

One thing you may notice about the [`BiodiversityObservationNetwork`](@ref) generated using [`SimpleRandom`](@ref) is that many of the points are clumped together. Many of the sampling algorithms in BiodiversityObservationNetworks aim to produce points that are _spatially balanced_, meaning they are well spread out across space, with little clumping. 

One such sampler is [`BalancedAcceptance`](@ref). Let's similarly make a [`BiodiversityObservationNetwork`](@ref) with 300 nodes that are spatially balanced.

```@example 1
bon = sample(BalancedAcceptance(150))
bonplot(bon, axistype=GeoAxis)
```

Much better! However, we are still missing some crucial things here. For example, what if we only want to select sites on land? This brings us to applying sampling algorithms to different _geometries_.

## Geometries in BiodiversityObservationNetworks.jl

SDMLayer, Polygon, Vectors of each, BONs themselves... the possibilities are increbidle.