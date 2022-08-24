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

    gradient()

end

dims, nl = (50, 50), 5
layers = [rand(MidpointDisplacement(), dims) for i in 1:nl]

optimize(layers, x->entropy(x))


function _squish(layers::Array{T, 3}, W::Matrix{T}) where {T <: AbstractFloat}
    return mapslices(x -> x * W, layers; dims = (2, 3))
end

function _squish(layers::Array{T, 3}, α::Vector{T}) where T <: AbstractFloat
    return reshape(mapslices(x -> x * α, layers; dims = (2, 3)), size(layers)[1:2]...)
end


heatmap(sl)