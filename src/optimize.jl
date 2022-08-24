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

function optimize(layers, loss; numtargets = 3, fixed_W = false, numsteps=10)
    numlayers = size(layers, 3)

    W = rand(numlayers, numtargets)
    α = rand(numtargets)

    η = 10^-4

    losses = zeros(numsteps)


    @showprogress for step in 1:numsteps
        grad = gradient(loss, W, α)
        ∂W, ∂α = η .* gradient(loss, W, α)
        W += ∂W
        α += ∂α

        losses[step] = loss(W,α)
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





model = (W, α) -> StatsBase.entropy(_squish(_squish(layers, W), α))/prod(size(layers[:,:,1]))
model(W, α)

optimize(layers, model; numsteps=10^4)


# Test run
x = _squish(layers, W)
typeof(x)
_squish(_squish(layers, W), α)


a = 10^-4 .* gradient(model, W, α)

numsteps = 10^5
losses = zeros(numsteps)
η = 10^-4
@showprogress for step in 1:numsteps
    grad = gradient(model, W, α)
    ∂W, ∂α = η .* gradient(model, W, α)
    W += ∂W
    α += ∂α

    losses[step] = model(W,α)


end

plot(1:length(losses), losses)

heatmap(model(W, α))