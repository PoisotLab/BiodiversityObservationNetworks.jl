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
using ProgressMeter

function _squish(layers::Array{T, 3}, W::Matrix{T}) where {T <: AbstractFloat}
    return convert(Array, slicemap(x -> x * W, layers; dims = (2, 3)))
end

function _squish(layers::Array{T, 3}, α::Vector{T}) where {T <: AbstractFloat}
    return slicemap(x -> x * reshape(α, (length(α), 1)), layers; dims = (2, 3))[:, :, 1]
end

squish(layers, W, α) = _squish(_squish(layers, W), α)


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

optimize(layers, model; numsteps = 10^4)
