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

function optimize(layers, simulator; numtargets = 3, fixed_W = false)
    numlayers = length(layers)

    W = rand(numlayers, numtargets)
    α = rand(numtargets)

    score = _squish(_squish(layers, Matrix(1.0I, 5, 3)), [0.3, 0.4, 0.3])
    loss = simulator(score)
    @info loss

    return gradient()
end

# ...?

using Statistics, StatsBase
using NeutralLandscapes
using Zygote, SliceMap

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

model = (W, α) -> StatsBase.entropy(_squish(_squish(layers, W), α))

# Test run
x = _squish(layers, W)
typeof(x)
_squish(_squish(layers, W), α)

model(W, α)

gradient(model, W, α)

heatmap(model(W, α))