"""
    What should the optimization API look like?


    We have a bunch of rasters that are layers that we 
    can to combine in a layer 'stack'.

    We have a weights matrix W which has `r` rows and
    `t` columns, where `r` is the number of layers in the stack,
    and `t` is the number of optimization targets.

    We have a vector α of length `t` which sums to 1.

    We want a function `optimize` which takes 
        (a) a combined Seeder/Refiner
        (b) an initial W and α
        (c) hyperparameters for optimization
        (d) a loss function comparing the sampled outcome to the 'true' state
    
    and uses Zygote's AD to optimize W and α
    to reduce a loss function that describes the 
    difference between the "true" metaweb and 
    sampled one. 

    In our context, the two targets are interaction classification
    and network topology, so we want a loss function that combines
    measures of these elements.
"""

struct Weights{F <: AbstractFloat}
    W::Matrix{F}
    α::Vector{F}
end

function optimize(layers, loss; targets = 3, learningrate = 1e-4, numsteps = 10)
    W = rand(size(layers, 3), targets)
    α = rand(targets)

    losses = zeros(numsteps)

    @showprogress for step in 1:numsteps
        ∂W, ∂α = learningrate .* gradient(loss, W, α)
        W += ∂W
        α += ∂α
        W = clamp.(W, 0, 1)
        α = clamp.(α, 0, 1)
        α ./= sum(α)

        losses[step] = loss(W, α)
    end
    return losses
end

# ...?

using Statistics, StatsBase
using ProgressMeter
using NeutralLandscapes
using Zygote, SliceMap
using ProgressMeter

function _squish(layers::Array{T, 3}, W::Matrix{T}) where {T <: AbstractFloat}
    return convert(Array, slicemap(x -> x * W, layers; dims = (2, 3)))
end

function _squish(layers::Array{T, 3}, α::Vector{T}) where {T <: AbstractFloat}
    return slicemap(x -> x * reshape(α, (length(α), 1)), layers; dims = (2, 3))[:, :, 1]
end

dims, nl, nt = (50, 50), 5, 3
W = rand(nl, nt)
α = rand(nt)
layers = zeros(dims..., nl)
for i in 1:nl
    layers[:, :, i] = rand(MidpointDisplacement(), dims)
end

model =
    (W, α) ->
        StatsBase.entropy(_squish(_squish(layers, W), α)) / prod(size(layers[:, :, 1]))

optimize(layers, model; numsteps = 10^4)