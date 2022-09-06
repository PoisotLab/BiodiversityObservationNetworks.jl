# # Selecting environmentally unique locations

# For some applications, we want to sample a set of locations with 
# most unique set of environmental coviarates. Another way to rephrase
# this problem is to say we want to find the set of points with the 
# _least_ covariance in their environmental values. 

using BiodiversityObservationNetworks
ENV["SDMLAYERS_PATH"] = "/home/michael/Data/RasterData/"
using SimpleSDMLayers
using NeutralLandscapes
using Plots


bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7)
temp, precip, seasonality, elevation = 
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 7; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 12; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 4; bbox...)),
    convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; bbox...))


function stack(layers::Vector{T}) where {T<:SimpleSDMLayer}
    # asset all layers are the same 
    mat = zeros(size(first(layers))..., length(layers))

    for (l,layer) in enumerate(layers)
        mat[:,:,l] .= broadcast(x->isnothing(x) ? NaN : x, layer.grid)

    end
    mat
end

layers = stack([temp,precip,seasonality,elevation])


# TODO this should be a extension of NL.jl,
# dispatch on the case that the kwarg `mask`
# is of type SDMLayer.
qcmask = fill(true, size(temp))                  
qcmask[findall(isnothing, temp.grid)] .= false   
uncert = rand(MidpointDisplacement(0.8), size(temp), mask=qcmask);
heatmap(uncert, aspectratio=1, frame=:box) 

candpts, uncert = uncert |> seed(BalancedAcceptance(α=0.01)); 
finalpts, uncert = refine(candpts, Uniqueness(;layers=layers), uncert)

heatmap(temp.grid)
scatter!([p[2] for p in candpts], [p[1] for p in candpts], c=:white)
scatter!([p[2] for p in finalpts], [p[1] for p in finalpts], c=:dodgerblue)