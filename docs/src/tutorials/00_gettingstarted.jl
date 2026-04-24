# # Getting Started with BiodiversityObservationNetworks.jl

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

# The core function is [`sample`](@ref). At minimum it needs a sampler and a
# domain. The simplest domain is a plain Julia `Matrix`:

mat = rand(50, 50)   # 50×50 grid of random values (imagine an elevation raster)

# and the simplest sampler is [`SimpleRandom`](@ref), which randomly selects sites without replacement.

result = sample(SimpleRandom(10), mat)

# [`sample`](@ref) returns a [`BiodiversityObservationNetwork`](@ref). It carries the selected
# site indices, their coordinates, and the values of auxiliary variables at the selected sites.


# TODO plotting with scatter (requires simple Makie extension)

scatter(result.coordinates; axis=(;aspect=1))


# When the domain is a single matrix (like `mat`), the auxiliary variables in the [`BiodiversityObservationNetwork`](@ref) 
# are simply the values of the original matrix at each selected.

# We can verify this by taking a look at the internals of the [`BiodiversityObservationNetwork`](@ref). 
# The sites are stored in a field called `sites`

result.sites

# and the auxiliary data associated which each site is stored in a matrix called `features`

result.features

# We can see that the first feature is the value of the domain at the first site

mat[result.sites[begin]]

# [`BiodiversityObservationNetwork`](@ref) also have a field called `coordinates`

result.coordinates

# At first, this may seem to be redundant as the same information is stored in `sites`,
# but this allows for storing both the Cartesian indices of a raster, and their corresponding 
# geospatial coordinates when using domains from [`SpeciesDistributionToolkit.jl`](). You can read
# more about the different types of domains [`here`](todo). 


# ## Masking sites

# Use the `mask` keyword to restrict sampling to a subset of the grid.
# A `true` value in the mask means the cell is *valid* (can be sampled):

mask = falses(50, 50)
mask[10:40, 10:40] .= true   # only sample the central region
result_masked = sample(SimpleRandom(10), mat; mask)

# We can plot the valid region in white and the invalid region in grey to see all 
# the coordinates fall in the valid region of the mask. 

heatmap(mask, colormap=[:grey50, :grey98])
scatter!(result_masked.coordinates)
current_figure()