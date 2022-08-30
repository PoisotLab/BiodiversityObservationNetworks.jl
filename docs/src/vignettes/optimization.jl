# # Optimization of layer weighting

# At the moment, we can only consider determinsitic
# loss functions because the other type of optimization 
# is hard and only sometimes converges ¯\_(ツ)_/¯
#
# For example, what if you want to optimize your
# spatial locations such that they have the least
# amount of covariance in environmental conditons-
# that is, they are as environmentally _unique_ as 
# possible. That is what we'll be building in this
# vignette.


using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots
using SimpleSDMLayers
using Zygote, SliceMap
using Statistics, StatsBase


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

#out = optimize(layers, model; numsteps = 10^4)
#heatmap(out)