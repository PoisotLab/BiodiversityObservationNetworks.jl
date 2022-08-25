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
        # Do some simulation with the current W and α
        candidatepts = _squish(_squish(layers, W), α) |> seed(BalancedAcceptance()) |> first
    
        truemetaweb, observedmetaweb = simulate_sampling(candidatepts)
        L = metaweb_loss(truemetaweb, observedmetaweb)

        ∂W, ∂α = learningrate .* gradient(loss, W, α)
        W += ∂W
        α += ∂α
        W = clamp.(W, 0, 1)
        α = clamp.(α, 0, 1)
        α ./= sum(α)


        numcolumns = size(W,2)
        for i in 1:numcolumns
            W[:,i] ./= sum(W[:,i])
        end

        losses[step] = loss(W, α)
    end
    return losses
end


using Statistics, StatsBase
using ProgressMeter
using NeutralLandscapes
using Zygote, SliceMap
using ProgressMeter
using Optim

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
        -StatsBase.entropy(_squish(_squish(layers, W), α)) / prod(size(layers[:, :, 1]))


function J(trueweb, obsweb)
    M, S = adjacency.([trueweb,obsweb])

    tp,tn,fp,fn = 0,0,0,0

    for i in eachindex(M)
        tp += M[i] == 1 && S[i] == 1
        tn += M[i] == 0 && S[i] == 0
        fp += M[i] == 0 && S[i] == 1
        fn += M[i] == 1 && S[i] == 0
    end
    return tp/(tp+fn) + tn/(tn+fp) - 1
end

function metaweb_loss(trueweb, obsweb) 
    # Compute some validation stats.
    # Topology using β-div from EN.jl
    # Interaction classification using Youden's J.

    #topological_loss = KGL01(βos(truemat, obsmat))
    interaction_loss = J(trueweb, obsweb)

    #@info "Topology Loss: $topological_loss"
    #@info "Interaction Loss: $interaction_loss"
    return interaction_loss #+ topological_loss
end 

metaweb_loss(truemat, obsmat)

losses = optimize(layers, model; numsteps = 10^4)

plot(1:length(losses), losses)


targets = 3 
W = rand(size(layers, 3), targets)
α = rand(targets)
layers

function stochastic_metaweb_loss(θ, layers; nlayers=3, ntargets=2)
    w_endpoint = nlayers*ntargets
    W, α = reshape(θ[1:w_endpoint], nlayers, ntargets), θ[w_endpoint+1:end]
    candidatepts = _squish(_squish(layers, W), α) |> seed(BalancedAcceptance()) |> first
    
    truemetaweb, observedmetaweb = simulate_sampling(candidatepts)
    L = metaweb_loss(truemetaweb, observedmetaweb)
    return L
end 

initθ = [vec(rand(5,3))..., vec(rand(3))...]

@time optimize(θ->stochastic_metaweb_loss(θ,layers; nlayers=nl, ntargets=nt), 
   initθ,
   ParticleSwarm(),
   Optim.Options(time_limit = 15.0))
)
