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


function optimize(layers, loss; targets = 3, learningrate = 1e-4, numsteps = 10)
    W = rand(size(layers, 3), targets)
    α = rand(targets)

    losses = zeros(numsteps)

    @showprogress for step in 1:numsteps
        ∂W, ∂α = learningrate .* gradient(loss, layers, W, α)
        W += ∂W
        α += ∂α
        W = clamp.(W, 0, 1)
        α = clamp.(α, 0, 1)
        α ./= sum(α)

        losses[step] = loss(W, α)
    end
    return losses
end
