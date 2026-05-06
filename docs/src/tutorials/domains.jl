# # Geospatial Domains in BONs.jl

# In most contexts, selecting spatial sampling sites requires choosing a set of locations from an explicitly geospatial domain, associated with geographical coordinates.
# BONs.jl supports both raster and vector geospatial domains through integration with [SpeciesDistributionToolkit.jl](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/) (SDT.jl).

# Here, we will explore how this integration works through examples.

# SDT.jl integration is included as a Julia _extension_, meaning each package must be installed and  loaded in order to enable functionality.

using BiodiversityObservationNetworks
using SpeciesDistributionToolkit

import Random #hide
Random.seed!(1234567890); #hide
using CairoMakie #hide
CairoMakie.activate!(; px_per_unit = 3) #hide

# # A Geospatial Raster Domain

# A natural extension of the `Matrix` domain, the simplest domain supported by BONs.jl, is using a geospatial raster. A raster is also just a matrix of values, but also containing additional metadata that associates each element of the matrix with a geographic location.

# BONs.jl supports this through the `SDMLayer` type from [SpeciesDistributionToolkit.jl](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/).

# Let's load an example SDMLayer to understand the basics of how the `SDMLayer` type works. Full documentation of what data sources are available in SDT, and how they can be manipulated, is available in the SDT.jl [documentation](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/).
# The below lines load a raster that contains the average annual temperature for Corsica. 

spatial_extent = (left = 8.412, bottom = 41.325, right = 9.662, top = 43.060) # the longitude/latitude extent of Corsica
temp = SDMLayer(RasterData(CHELSA1, AverageTemperature); spatial_extent...,)

# We can plot it using `heatmap` from Makie to visualize


f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, temp, colormap=:OrRd)
current_figure() #hide


# Sampling from this `SDMLayer` works exactly like sampling from a matrix. For example, we can use the [`SimpleRandom`](@ref) as follows:

bon = sample(SimpleRandom(), temp)

# Which we can then plot 

f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, temp, colormap=:OrRd)
scatter!(ax, bon, color=:white, strokewidth=1, strokecolor=:black)
current_figure() #hide

# Note that there are regions in the raster (the water of surrouding Corsica) that have no value. Samplers will automatically avoid sampling sites in those regions --- only pixels with valid data are considered for sampling. 

# ## Masking Geospatial Rasters

# As we've just seen, regions without valid pixel data are automatically ignored. However, we can mask addition pixels using the `mask` keyword.

# For a SDMLayer domain, the type of object we pass to `mask` can either be *(1)* a matrix with the same size as the `SDMLayer` domain. or *(2)* another `SDMLayer` with the *same* size, extent, and resolution as the `SDMLayer` used as the domain, 

# Here we will provide an example of both forms to mask out the center region of Corsica.

# ### Masking a `SDMLayer` with a matrix 

# Let's construct our mask by first creating a matrix of all ones that is the same size as our `SDMLayer`.

matrix_mask = ones(size(temp));

# And now set a box in the middle to 0s to make it "off-limits" for the sampler

matrix_mask[100:150, 50:125] .= false;

# Now we can sample 

masked_bon = sample(SimpleRandom(), temp, mask = matrix_mask)

# And plot to verify points are all outside the masked region.

temp2 = copy(temp) #hide 
temp2.indices[findall(iszero, matrix_mask)] .= 0 #hide

f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, temp2, colormap=:OrRd)
scatter!(ax, masked_bon, color=:white, strokewidth=1, strokecolor=:black)
current_figure() #hide

# ::: note Hiding masked regions in the SDMLayer
#
# Note that we have done a bit of trickery on the SDMLayer to hide the masked pixels in the plot.
# This was done by running 
# ```julia
# temp2 = copy(temp) 
# temp2.indices[findall(iszero, matrix_mask)] .= 0 
# ```
# and plotting `temp2`, which we sneakily hid from the code above to avoid confusion and unnecessary detail.
#
# Note this has *no effect* on the actual sampling of the sites, it is simply for visualization purposes to verify no sampled sites fall within the mask.
# 
# :::

# ### Masking a `SDMLayer` with another `SDMLayer`

# `SDMLayer`s can also be used as masks, assuming they have the *same size*, *extent*, and *crs* as the domain.

# *SDMLayer masks only mask out the regions that are set as invalid pixels, regardless of the pixel value*.
# Unlike a matrix mask, the whether the value of an SDMLayer's pixel is *not what is used to determine if a site is masked*.
# `SDMLayer`s have a separate field called `indices`, which determine what pixels are valid. *The `indices` field is what is used to mask sampling sites*

# To construct an `SDMLayer` mask, we start by constructing a `copy` of our temp

layer_mask = copy(temp)

# The validity of each pixel is stored in the `indices` field of an SDM. To construct a mask, we can run the lines

layer_mask.indices[100:150, 50:125] .= 0; # set center of Corsica to 0

# which will remove the center region of Corsica from consideration.

# We can visualize the new layer mask

f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, layer_mask, colormap=:OrRd)
current_figure() #hide

# We can then pass the `layer_mask` to the `mask` keyword argument as normal:

masked_bon = sample(SimpleRandom(), temp, mask = layer_mask)

# and plot to verify the masking works


f = Figure()
ax = Axis(f[1,1], aspect=DataAspect())
hm = heatmap!(ax, layer_mask, colormap=:OrRd)
scatter!(ax, masked_bon, color=:white, strokewidth=1, strokecolor=:black)
current_figure() #hide


# ## Custom Inclusion Probabilities with Geospatial Rasters

# When using custom inclusion probabilities with an `SDMLayer` domain, things work very similarly to masking.
# Either (1) a matrix of the same size or (2) an `SDMLayer` with the same size, extent and crs can be used as inclusion probabilities. 

# ### A inclusion probabilities matrix with an `SDMLayer` domain.

# Let's construct inclusion probabilities where the right side is more likely to be included, as we did in the first tutorial.

inclusion_probabilies = [1.15^j for i in 1:size(temp,1), j in 1:size(temp, 2)];

# ::: note x/y vs. longitude/latitude
#
# You may notice above that to make inclusion increase as we go right, we use `1.15^j`, where as in the first tutorial we used `1.15^i`.
#  
# This is because in a matrix, `i` (the first axis) typically corresponds to the horizontal axis, but in a raster, `j` (the second axis) typically corresponds to the horizontal/longitude axis. 
# 
# This is a result of decisions that people made before many of us were born and that we will all likely have to live with until we are all dead. That's just how it is.
#
# :::

# We can then use this with the `inclusion` keyword argument

inclusion_bon = sample(SimpleRandom(), temp, inclusion = inclusion_probabilies)

# and visualize to confirm it worked


f = Figure()
ax = Axis(f[1,1])
heatmap!(ax, temp, colormap=[:grey80])
scatter!(inclusion_bon, color=:dodgerblue)
current_figure() #hide
 

# # Using a geospatial vector as a domain

# SpeciesDistributionToolkit also supports various types of geospatial vector data (e.g. polygons) as both domains and rasters

# it rasterizes it because all samplers work on discrete domains.
# the resolution it rasterizes it at is customizable.

# polygons can also be used as masks for SDMLayers, or other polygons.

# if you want to use custom inclusion probabilities, it's possible using a polygon directly.
# its going to be more straightforward to construct an SDMLayer masked by the polygon with mask! 
# and base the inclusion probabilities on that. 





