# # Selecting environmentally unique locations

# For some applications, we want to sample a set of locations that cover a broad
# range of values in environment space. Another way to rephrase this problem is
# to say we want to find the set of points with the _least_ covariance in their
# environmental values. 

# To do this, we use a `BONRefiner` called `Uniqueness`. We'll start by loading
# the required packages. 

using BiodiversityObservationNetworks
using SimpleSDMLayers
using StatsBase
using NeutralLandscapes
using Plots

# !!! warning "Consider setting your SDMLAYERS_PATH" When accessing data using
#     `SimpleSDMLayers.jl`, it is best to set the `SDM_LAYERSPATH` environmental
#     variable to tell `SimpleSDMLayers.jl` where to download data. This can be
#     done by setting `ENV["SDMLAYERS_PATH"] = "/home/user/Data/"` or similar in
#     the `~/.julia/etc/julia/startup.jl` file. (Note this will be different
#     depending on where `julia` is installed.)

bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7);
temp, precip, seasonality, elevation = 
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 7; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 12; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 4; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; bbox...));

# Now we'll use the `stack` function to combine our four environmental layers
# into a single, 3-dimensional array, which we'll pass to our `Uniqueness` refiner.

layers = stack([temp,precip,seasonality,elevation]);


# this requires NeutralLandscapes v0.1.2
uncert = rand(MidpointDisplacement(0.8), size(temp), mask=temp);
heatmap(uncert, aspectratio=1, frame=:box) 

# Now we'll get a set of candidate points from a BalancedAcceptance seeder that
# has no bias toward higher uncertainty values.
candpts, uncert = uncert |> seed(BalancedAcceptance(numpoints=100, α=0.0)); 

# Now we'll `refine` our `100` candidate points down to the 30 most
# environmentally unique.
finalpts, uncert = refine(candpts, Uniqueness(;numpoints=30, layers=layers), uncert)

heatmap(uncert)
scatter!([p[2] for p in candpts], [p[1] for p in candpts], fa=0.0, msc=:white, label="Candidate Points")
scatter!([p[2] for p in finalpts], [p[1] for p in finalpts], c=:dodgerblue, msc=:white, label="Selected Points")