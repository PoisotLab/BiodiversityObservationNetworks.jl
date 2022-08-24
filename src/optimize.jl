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
"""

struct Weights{F<:AbstractFloat}
    W::Matrix{F}
    α::Vector{F}
end

function optimize(layers, proposer, loss)

end
