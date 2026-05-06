# # Getting Started with BiodiversityObservationNetworks.jl

# BiodiversityObservationNetworks.jl (BONs.jl) is a Julia package for selecting sites for spatial sampling, typically in the context of ecological or biodiversity sampling design. 
# The purpose of this package is to provide a high-level, extensible, modular interface to the selection of sampling point for biodiversity processes in space.
# BONs.jl implements a variety of algorithms for site selection, including stratified and spatially balanced sampling, and enables adaptive sampling based on model-based estimates of uncertainty or the locations of legacy sampling sites. 
# It also includes a variety of utilities for quantifying a sample's spatial balance and how representative a sample is of auxiliary environmental variables.

using BiodiversityObservationNetworks
import Random
using CairoMakie
Random.seed!(1234567890); #hide
CairoMakie.activate!(; px_per_unit = 3) #hide

# In this tutorial, we will cover the basics of how to use BiodiversityObservationNetworks.jl. We'll start by loading the package.

# ## Installing and loading

# Install the package from the Julia REPL:
# ```julia
# using Pkg; Pkg.add("BiodiversityObservationNetworks")
# ```
# Then load it:
# ```julia
# using BiodiversityObservationNetworks
# ```

# ## Your first sample

# The core function in BONs.jl is [`sample`](@ref). [`sample`](@ref) works by taking at minimum a sampling algorithm (a type of [`BONSampler`](@ref)) and a *domain*, which represents the region from which to select sampling points. 

# The core pattern is to call `sample(algorithm, domain)`, with additional arguments available depending on the specific algorithm used.  

# BONs.jl supports a variety of domains, including raster and vector data. The simplest domain is a plain Julia `Matrix`:

mat = rand(50, 50);   # 50×50 grid of random values (e.g. an elevation raster with very weird topography)

# The simplest sampler is [`SimpleRandom`](@ref), which randomly selects sites without replacement.

result = sample(SimpleRandom(), mat)

# [`sample`](@ref) returns a [`BiodiversityObservationNetwork`](@ref), which stores the selected
# site indices, their coordinates, and (if available) the values of auxiliary variables at the selected sites.

# By default, all algorihtms will select 50 sites. This is changed simply by passing an integer value to the sampler, e.g.

result = sample(SimpleRandom(10), mat)

# will yield 10 sites. 

# ## Visualizing Results

# To plot the results, we use Makie, a very feature-rich package in Julia for data visualization.

# BONs.jl functionality for Makie relies on a Julia extension, meaning it is only activated if Makie is loaded in the same environment.

# Makie can be installed and loaded similar to above. Makie provides various "backends", for different types of plotting. Here we will used CairoMakie, which is the backend for publication-quality graphics. You can learn more about Makie and the different available backends [here](https://makie.org/website/)

# ```julia
# using CairoMakie
# ```

# A BON can then be visualized with the `scatter` method. 

scatter(result)

# All typical keyword arguments can be passed to `scatter` to customize the result, e.g.

scatter(result; color = :green, marker = :star5, markersize = 15, axis=(;aspect=1))

# ## Understanding the BiodiversityObservationNetwork type

# The [`BiodiversityObservationNetwork`](@ref) type stores a variety of information about the sample. 

# The first is `sites`, which are the indices in the original domain for the selected sites (typically as a vector of `CartesianIndex`s)

result.sites

# [`BiodiversityObservationNetwork`](@ref) also have a field called `coordinates`

result.coordinates

# At first, this may seem to be redundant as the same information is stored in `sites`,
# but this allows for storing both the Cartesian indices of selected sites a raster, and their corresponding 
# geospatial coordinates when using supported geospatial domains from [`SpeciesDistributionToolkit.jl`](TODO). You can read
# more about the different types of domains [`here`](TODO). 

# When the domain is a single matrix (like `mat`), the auxiliary variables in the [`BiodiversityObservationNetwork`](@ref) 
# are simply the values of the original matrix at each selected, which are stored in a matrix called `features`

result.features

# We can see that the first feature is the value of the domain at the first site

mat[result.sites[begin]]

# When using multiple rasters/matrices as a domain, the auxiliary variables are the values across each of those rasters.
# For example, 

multilayer_domain = [rand(30,20) for i in 1:5]
multilayer_result = sample(SimpleRandom(10), multilayer_domain)

# Will yield a matrix of auxiliary variables of dimension 5 x 10 (for 5 auxiliary variables at each of 10 sites) 

multilayer_result.features

# ## Masking sites

# Use the `mask` keyword to restrict sampling to a subset of the grid.
# A `true` value in the mask means the cell is *valid* (i.e. can be included in a sampled result):

# As an example, lets build a mask that restricts points to the center 30x30 region of a 50x50 matrix.

mask = falses(50, 50)
mask[10:40, 10:40] .= true   # only sample the central region
result_masked = sample(SimpleRandom(10), mat; mask)

# We can plot the valid region in white and the invalid region in grey to verify that all 
# the coordinates fall in the valid region of the mask. 

heatmap(mask, colormap=[:grey50, :grey98])
scatter!(result_masked.coordinates)
current_figure()

# ## Custom Inclusion Probabilities 

# Up until this point, we have been using the [`SimpleRandom`](@ref) sampler, which draws sampling locations with equal probability, without replacement.
# For a variety of reasons, we may want the initial probability that a site is included to vary, so some locations are more likely to be included than others. 
# Many sampling algorithms, including [`SimpleRandom`](@ref), support this functionality. 

# For example, lets make it so the inclusion probability increases as we move from left to right across the domain. This can be done via

inclusion_probability = [1.1^i for i in 1:50, j in 1:50]

# Let's plot this matrix to verify this is what we get

heatmap(inclusion_probability)

# Now we can sample with these custom inclusion probabilities using the `inclusion` keyword argument

result = sample(SimpleRandom(), inclusion_probability, inclusion = inclusion_probability)

# and very can verify the selected points are skewed more toward the right side of the domain

scatter(result)

# Note that the inclusion probability matrix doesn't _have_ to be the domain. 
# Different inclusion probabilities and domains can be used as long as they are compatible,
# meaning they are equally sized matrices, or a SDMLayer with matching size, extent, and crs if using the SDMLayers extension.

# For example

result = sample(SimpleRandom(), rand(50, 50); inclusion = inclusion)

# is also valid.

# ## Next Steps

# This covers the basic functionality of BONs.jl. More advanced functionality is explored in the following tutorials, which we recommend following in the below order:
#  - [Geospatial Domains with SpeciesDistributionToolkit](./tutorials/domains)
#  - [Multistage Samplers](./tutorials/multistage)
#  - [Targeting unique climates or regions with high climate velocity](./tutorials/climate) 
#  - [Evaluating selected sites](./tutorials/evaluation)
#  - [Including Legacy Sampling Sites](./tutorials/legacy)
#  - [Adaptive Sampling of Species Distributions](./tutorials/adaptive)

# In addition, full descriptions of each of the *supported algorithms for sampling* can be found [here](TODO).

# Description of the various utilities in BONs.jl can be found [here](), and a design document describing how the internals of the package are designed (primarily aimed for contributors to the package) can be found [here](TODO).
