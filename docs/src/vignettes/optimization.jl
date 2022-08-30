# # Optimization of layer weighting

# At the moment, we can only consider determinsitic
# loss functions because the other type of optimization 
# is hard and only sometimes converges ¯\_(ツ)_/¯
#
# For example, what if you want to optimize your
# spatial locations such that they have the least
# amount of covariance in environmental conditons-
# that is, they are as environmentally _unique_ as 
# possible. That is what we'll be building in this
# vignette.


using BiodiversityObservationNetworks
using NeutralLandscapes
using Plots
using SimpleSDMLayers
using Optim
using Zygote, SliceMap
using Statistics, StatsBase

using Optim


dims, nl, nt = (50, 50), 5, 3
W = rand(nl, nt)
α = rand(nt)
layers = zeros(dims..., nl)
for i in 1:nl
    layers[:, :, i] = rand(MidpointDisplacement(), dims)
end

function covariance_map(layers)
    #candidatepts = squish(layers, W, α) |> seed(BalancedAcceptance()) |> first

    covarmap = zeros(size(layers,1),size(layers,2))
    ind = CartesianIndices((1:size(covarmap,1),1:size(covarmap,2)))
    for i in ind
        vi = layers[i[1], i[2],:]
        for j in ind
            if i != j
                vj =  layers[j[1],j[2],:]
                covarmap[i] += abs(cov(vi,vj))
            end
        end
    end
    return covarmap
end 

initθ = [vec(rand(5,3))..., vec(rand(3))...]


heatmap(covariance_map(layers))


out = @time Optim.optimize(θ->covariance_loss(θ,layers; nlayers=nl, ntargets=nt), 
   initθ,
   ParticleSwarm(),
   Optim.Options(time_limit = 15.0))




#out = optimize(layers, model; numsteps = 10^4)
#heatmap(out)