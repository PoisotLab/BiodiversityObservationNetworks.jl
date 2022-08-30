# # Optimization of layer weighting

# At the moment, we can consider determinsitic
# loss functions.

using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots


# ...?

using Statistics, StatsBase
using ProgressMeter
using NeutralLandscapes
using Zygote, SliceMap


dims, nl, nt = (50, 50), 5, 3
W = rand(nl, nt)
α = rand(nt)
layers = zeros(dims..., nl)
for i in 1:nl
    layers[:, :, i] = rand(MidpointDisplacement(), dims)
end

model =
    (layers, W, α) ->
        StatsBase.entropy(squish(layers,W,α)) / prod(size(layers[:, :, 1]))

out = optimize(layers, model; numsteps = 10^4)
heatmap(out)