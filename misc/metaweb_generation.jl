using Distributions: BetaBinomial
using EcologicalNetworks: nichemodel
using BiodiversityObservation 

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
    richness=50,
    numobs = 100)

    metaweb = flexiblelinks(richness)
    ra = getabundances(metaweb)

    for coord in coords
        observe(metaweb, coords, ra, numobs)
    end
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

function observe(metaweb, coord, relativeabundances, numobservations, )
    counts = zeros(length(relativeabundances))
    for i in 1:numobservations
        ind = rand(Categorical(relativeabundances))
        counts[ind] += 1
    end

    observed_metaweb = zeros(size(adjacency(metaweb)))

    S = richness(metaweb)

    for i in 1:S, j in 1:S
        if metaweb[i,j] == 1 && observations[i] > 0 && observations[j] > 0
            observed_metaweb[i,j] = 1
        end
    end
    return UnipartiteNetwork(observed_metaweb)
end



simulate_sampling()