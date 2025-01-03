# Selecting environmentally unique locations

For some applications, we want to sample a set of locations that cover a broad
range of values in environment space. Another way to rephrase this problem is to
say we want to find the set of points with the _least_ covariance in their
environmental values.  

To do this, we use a `BONRefiner` called `Uniqueness`. We'll start by loading the required packages. 

```@example 1
using BiodiversityObservationNetworks
using SpeciesDistributionToolkit
using StatsBase
using NeutralLandscapes
using CairoMakie
```

!!! warning "Consider setting your SDMLAYERS_PATH" 
    When accessing data using `SimpleSDMDatasets.jl`, it is best to set the `SDM_LAYERSPATH` environmental variable to tell `SimpleSDMDatasets.jl` where to download data. This can be done by setting `ENV["SDMLAYERS_PATH"] = "/home/user/Data/"` or similar in the `~/.julia/etc/julia/startup.jl` file. (Note this will be different depending on where `julia` is installed.)

```@example 1
bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7);
temp, precip, elevation = 
    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, AverageTemperature); bbox...)),
    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, Precipitation); bbox...)),
    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, Elevation); bbox...));
```

Now we'll use the `stack` function to combine our four environmental layers into a single, 3-dimensional array, which we'll pass to our `Uniqueness` refiner.

```@example 1
layers = BiodiversityObservationNetworks.stack([temp,precip,elevation]);
```

```@example 1
uncert = rand(MidpointDisplacement(0.8), size(temp), mask=temp);
heatmap(uncert) 
```

Now we'll get a set of candidate points from a BalancedAcceptance seeder that has no bias toward higher uncertainty values.

```@example 1
candpts = seed(BalancedAcceptance(numsites=100)); 
```

Now we'll `refine` our `100` candidate points down to the 30 most environmentally unique.

```@example 1
finalpts = refine(candpts, Uniqueness(;numsites=30, layers=layers))
heatmap(uncert)
scatter!([p[1] for p in candpts], [p[2] for p in candpts], color=:white)
scatter!([p[1] for p in finalpts], [p[2] for p in finalpts], color=:dodgerblue, msc=:white)
current_figure()
```
