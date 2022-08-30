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

function _squish(layers::Array{T, 3}, W::Matrix{T}) where {T <: AbstractFloat}
    return convert(Array, slicemap(x -> x * W, layers; dims = (2, 3)))
end

function _squish(layers::Array{T, 3}, α::Vector{T}) where {T <: AbstractFloat}
    return slicemap(x -> x * reshape(α, (length(α), 1)), layers; dims = (2, 3))[:, :, 1]
end

"""
    squish(layers, W, α)

Takes a set of `n` layers and squishes them down
to a single layer.

For a coordinate in the raster (i,j), denote the vector of
values across all locations at that coordinate v⃗ᵢⱼ. The value
at that coordinate in squished layer, s⃗ᵢⱼ, is computed in two steps.

**(1):** First we apply a weights matrix, `W``, with `n` rows and `m` columns (`m` < `n`), to
reduce the initial `n` layers down to a set of `m` layers, each of which corresponds 
to a particular target of optimization. For example, we may want to propose sampling 
locations that are optimized to best sample  abalance multiple criteria, like (a) the 
current distribution of a species and (b) if that distribution is changing over time.

Each entry in the weights matrix `W` corresponds to the 'importance' of the layer
in the corresponding row to the successful measurement of the target of the corresponding
column. As such, each column of `W` must sum to 1.0.

For each location, the value of the condensed layer `tᵢ`, corresponding to target `i`, at 
coordinate (i,j) is given by the dot product of v⃗ᵢⱼ and the `i`-th column of `W`.

**(2):** Apply a weighted average across each target layer. To produce the final output layer,
we apply a weighted average to each target layer, where the weights are provided in the vector α⃗
of length `m`.

The final value of the squished layer at (i,j) is given by s⃗ᵢⱼ = ∑ₓ αₓ*tᵢⱼ(x), where tᵢⱼ(x) is
the value of the x-th target layer at (i,j).
"""
squish(layers, W, α) = _squish(_squish(layers, W), α)
