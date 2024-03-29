@info "Loading BONs.jl support for Zygote.jl ..."


function BiodiversityObservationNetworks.optimize(layers, loss; targets = 3, learningrate = 1e-4, numsteps = 10)
    ndims(layers) == 3 || throw(ArgumentError("Layers must be a 3-dimensional array"))
    typeof(loss) <: Function || throw(ArgumentError("`loss` must be a function"))
    learningrate > 0.0 || throw(ArgumentError("learningrate must be positive"))
    numsteps > 0 || throw(ArgumentError("numsteps must be positive"))

    W = rand(size(layers, 3), targets)
    for (i,c) in enumerate(eachcol(W))
        W[:,i] .= c ./ sum(c)
    end

    α = rand(targets)
    α ./= sum(α)
    
    losses = zeros(numsteps)

    @showprogress for step in 1:numsteps
        dL, ∂W, ∂α = learningrate .* Zygote.gradient(loss, layers, W, α)
        W += ∂W
        α += ∂α

        W = clamp.(W, 0, 1)
        for (i,c) in enumerate(eachcol(W))
            W[:,i] .= c ./ sum(c)
        end

        α = clamp.(α, 0, 1)
        α ./= sum(α)

        losses[step] = loss(layers, W, α)
    end
    return W,α,losses
end
