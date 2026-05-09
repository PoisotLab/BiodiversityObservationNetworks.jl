# Case Study: Adaptive Sampling of Species Distributions

using Pkg

logistic(x, α, β) = 1 / (1 + exp((x - β) / α))
logistic(α, β) = (x) -> logistic(x, α, β)

@kwdef struct LogisticNiche
    centers = [rand(Uniform(0.2,0.8)) for i in 1:19]
    shapes = [rand(Normal(0,0.4)) for i in 1:19]
end

function predict(niche::LogisticNiche, x)
    predictors = [logistic(x[i], niche.shapes[i], niche.centers[i]) for i in eachindex(x)]
    prod(predictors)
end