# ===============================================
#      
#   This file contains methods to simulate the 
#   process of observing a set of local ecological
#   networks at a site of proposed coordinates. It
#   then returns the "true" metaweb, and the observed
#   metaweb, which can then be fed into a loss function
#   that measures the difference between the two using 
#   a handfull of metrics. 
# 
#   What should these loss metrics be? Necessary conversation
#   to have. Just like |κ_true - κ_obs | and | connectance_true - connectance_obs|?
# 
#   
# 
# ===============================================

using Distributions
using EcologicalNetworks
using BiodiversityObservationNetworks

function flexiblelinks(S)
    # MAP estimates from Macdonald et al 2020
    ϕ = 24.3
    μ = 0.086

    n = S^2 - (S-1)
    α = μ*ϕ
    β = (1-μ)*ϕ
    
    L = 0
    while L < 1 || (L > 0.5*S^2)
        L = rand(BetaBinomial(n, α, β))
    end
    return nichemodel(S, L)
end


function simulate_sampling(
    coords::Vector{CartesianIndex};
    speciesrichness=50,
    numobs = 10,
    metaweb=flexiblelinks(speciesrichness),
    sz = (50,50)
)
    ra = getabundances(metaweb)
    occurrence_maps = []
    for i in 1:richness(metaweb)
        thismap = classify(rand(MidpointDisplacement(), sz), [0.25, 1.]) .== 1 
        push!(occurrence_maps, thismap)
    end 
    locals = [UnipartiteNetwork(rand(Bool,5,5)) for _ in 1:length(coords)]
    for (i,coord) in enumerate(coords)
        locals[i] = observe(metaweb, coord, ra, numobs, occurrence_maps)
    end

    observed_metaweb = reduce(∪, locals)
    return metaweb, observed_metaweb
end 


function getabundances(metaweb)
    S = richness(metaweb)
    # By trophic levels (possible but slow): 
    # Z = 2
    # abundances = Z.^[trophdict["s$i"]-1 for i in 1:S]
    # abundance_dist = abundances ./ sum(abundances)

    # By lognormal dist:
    abundances = rand(LogNormal(),S)
    return  abundances ./ sum(abundances)
end


function observe(metaweb, 
    coord, 
    relativeabundances, 
    numobservations, 
    occurrence_maps
)

    counts = zeros(length(relativeabundances))

    occurrence = [occurrence_maps[i][coord] == 1 for i in eachindex(relativeabundances)]

    realized_ra = similar(relativeabundances)
    for i in eachindex(relativeabundances)
        realized_ra[i] = occurrence[i] == 1 ? relativeabundances[i] : 0
    end

    observed_net = zeros(size(adjacency(metaweb)))

    sum(realized_ra) == 0 && return UnipartiteNetwork(Bool.(observed_net))

    realized_ra = [i/sum(realized_ra) for i in realized_ra]
    for i in 1:numobservations
        ind = rand(Categorical(realized_ra))
        counts[ind] += 1
    end


    S = richness(metaweb)

    for i in 1:S, j in 1:S
        if metaweb[i,j] == 1 && counts[i] > 0 && counts[j] > 0
            observed_net[i,j] = 1
        end
    end
    return UnipartiteNetwork(Bool.(observed_net))
end

metaweb = UnipartiteNetwork(
          Bool[0 0 0;
           0 0 0;
           1 1 0])



pts = rand(50,50) |> seed(BalancedAcceptance()) |> first
simulate_sampling(pts, metaweb=metaweb)


truemat, obsmat = simulate_sampling(pts)




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
    @info "Interaction Loss: $interaction_loss"
    return 1-interaction_loss #+ topological_loss
end 

metaweb_loss(truemat, obsmat)

losses = optimize(layers, model; numsteps = 10^4)

plot(1:length(losses), losses)


targets = 3 
W = rand(size(layers, 3), targets)
α = rand(targets)
layers

function stochastic_metaweb_loss(θ, layers; metaweb=flexiblelinks(20),nlayers=3, ntargets=2)
    @info θ
    w_endpoint = nlayers*ntargets
    W, α = reshape(θ[1:w_endpoint], nlayers, ntargets), θ[w_endpoint+1:end]
    candidatepts = _squish(_squish(layers, W), α) |> seed(BalancedAcceptance(α=3.)) |> first
    
    truemetaweb, observedmetaweb = simulate_sampling(candidatepts, metaweb=metaweb)
    L = metaweb_loss(truemetaweb, observedmetaweb)
    return L
end 

initθ = [vec(rand(5,3))..., vec(rand(3))...]

out = @time optimize(θ->stochastic_metaweb_loss(θ,layers; nlayers=nl, ntargets=nt), 
   initθ,
   ParticleSwarm(),
   Optim.Options())
)