using Distributions: BetaBinomial
using EcologicalNetworks: nichemodel

function flexiblelinksmodel(S)
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

