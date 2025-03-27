# Getting Started with BiodiversityObservationNetworks.jl

In this tutorial, we will cover the basics of how to use BiodiversityObservationNetworks.jl. We'll start by loading the package.

```@ansi 1
using BiodiversityObservationNetworks
```

## "Hello World" in BiodiversityObservationNetworks.jl

The primary use of BiodiversityObservationNetworks is for generating [`BiodiversityObservationNetwork`](@ref) using a variety of point-selection algorithms. Generating networks, regardless of the specific algorithm chosen, is done using the [`sample`](@ref) method. The simplest point-selection algorithm is [`SimpleRandom`](@ref), where each location in space (what is meant be location? More on this in the next section) has an equal probability of inclusion.

```@ansi 1
bon = sample(SimpleRandom())
```

By default, [`SimpleRandom`](@ref) (and every other sampler), chooses 50 points. Without any other arguments, [`sample`](@ref) chooses points from a raster than covers the entire globe, with each pixel representing a 1˚ by 1˚ region.

When a version of the [Makie](https://docs.makie.org/v0.22/) package for data visualization is loaded, we can use built-in functions to visualize the network. Let's use [`GeoMakie`](https://geo.makie.org/v0.7.9/ ), which specifically enables plotting of geographic plots within Makie, with the `CairoMakie` backend. (Learn mroe about Makie backends [here](https://docs.makie.org/stable/explanations/backends/backends#What-is-a-backend))


```@ansi 1
using CairoMakie, GeoMakie
```

Now that these packages are loaded, we can use the [`bonplot`](@ref) method to visualize the generated [`BiodiversityObservationNetwork`](@ref). We pass `GeoAxis` to the keyword argument `axistype` to ensure [`bonplot`](@ref) uses `GeoMakie`.


```@ansi 1
bonplot(bon, axistype = GeoAxis)
```


## 