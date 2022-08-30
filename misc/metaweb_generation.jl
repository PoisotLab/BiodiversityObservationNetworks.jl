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

