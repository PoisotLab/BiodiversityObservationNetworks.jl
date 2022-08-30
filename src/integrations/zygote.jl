@info "Loading BONs.jl support for Zygote.jl ..."


function BiodiversityObservationNetworks.optimize(layers, loss; targets = 3, learningrate = 1e-4, numsteps = 10)
    W = rand(size(layers, 3), targets)
    α = rand(targets)

    losses = zeros(numsteps)

    @showprogress for step in 1:numsteps
        ∂W, ∂α = learningrate .* Zygote.gradient(loss, layers, W, α)
        W += ∂W
        α += ∂α
        W = clamp.(W, 0, 1)
        α = clamp.(α, 0, 1)
        α ./= sum(α)

        losses[step] = loss(W, α)
    end
    return losses
end
