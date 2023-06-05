# kldivergence.jl
# -------------------------------------- 
#
# 
# The idea for this method is, given a stack of rasters arbitrary size with
# values that have mean μᵢ and variance σᵢ for each layer `i`, to find the `n`
# cells to sample that have mean (μₛᵢ) and variance (σ²ₛᵢ) for each layer, such
# that the mean Kullback–Leibler divergence between Normal(μᵢ, σᵢ²) and
# Normal(μₛ,σ²ₛ) is minimized across all `i`. 
#
# 
# There are two different ideas for how we can implement this: 
#
# 1. Using autodifferentiation: The points should be selected by applying a
#    single layer perceptron to the raster followed by softmax, and pick the top
#    `n` values as the candidate points. The weights and biases for the single
#    layer can be optimized by taking gradients using Zygote. 
#
# 2. Using stochastic optimization: The points can be initially proposed to be
#    balanced, and using simulated annealing to do the optimization. The
#    proposal function would be moving a single point at a time by an offset
#    with a kernel that has distance proportional to current temperature. 
#
#


using Distributions
using NeutralLandscapes
using BiodiversityObservationNetworks
using ProgressMeter

using Optim

abstract type CoolingSchedule end 
decay(s::CoolingSchedule) = s.decay

Base.@kwdef struct ExponentialCooling <: CoolingSchedule 
    decay = 0.0001
end 
(ec::ExponentialCooling)(T₀, n) = exp(-decay(ec)*(n-1)) * T₀

Base.@kwdef struct GeometricCooling <: CoolingSchedule 
    decay = 0.9999
end 
(gc::GeometricCooling)(T₀, n) = (decay(gc)^(n-1)) * T₀

Base.@kwdef struct InverseLogCooling <: CoolingSchedule
    decay = 1.
end
(ilc::InverseLogCooling)(T₀, n) =  T₀ / (1+ decay(ilc)*log(n))


Base.@kwdef mutable struct KLSimulatedAnnealing{I<:Integer, F<:AbstractFloat} <: BONRefiner
    numpoints::I = 50
    layers::Array{F, 3} = rand(100,100,3)
    cooling::CoolingSchedule = ExponentialCooling()
    T₀::F = 10.
    proposal_kernel = Exponential(5)
    steps::I = 1000
    multivariate::Bool = false
end 

temperature_schedule(klsa::KLSimulatedAnnealing{I,F}) where {I,F} =  
[klsa.cooling(klsa.T₀,i) for i in 1:klsa.steps]

function move_probability(KL_proposed, KL_now,T)
    return exp(-(KL_proposed - KL_now) / T)
end 

function _clip_inbounds(coord, sz)
    x,y = sz 
    coord[1] <= 0 && coord[2] <= 0 && return CartesianIndex(1,1)
    coord[1] > x && coord[2] > y && return CartesianIndex(x,y)
    coord[1] > x && coord[2] <= 0 && return CartesianIndex(x,1)
    coord[2] > y && coord[1] <= 0 && return CartesianIndex(1,y)
    coord[1] <= 0 && return CartesianIndex(1, coord[2])
    coord[2] <= 0 && return CartesianIndex(coord[1], 1)
    coord[2] > y && return CartesianIndex(coord[1], y)
    coord[1] > x && return CartesianIndex(x, coord[2])
    return coord
end

function propose(klsa::KLSimulatedAnnealing, coord)
    r = rand(klsa.proposal_kernel)
    θ = rand(Uniform(0, 2π))
    Δx, Δy = Int32.(floor.([r*cos(θ), r*sin(θ)]))
    _clip_inbounds(coord + CartesianIndex(Δx, Δy), size(klsa.layers)[1:2])
end

# takes two vectors of CartesianIndices and mutates one of them.
# have to explore ways more args can be passed via SA 
function propose_mutating!(x_proposed, x_current)
    
end

function jensen_shannon_distance(P,Q; kwargs...)
    # there is also occasional numerical instability with very divergent distributions
    M = MixtureModel([P,Q], [0.5,0.5])
    div = 0.5kldivergence(P,M; kwargs...) + 0.5kldivergence(Q,M;kwargs...)
    sqrt(div / log(2))
end

function _multivariate_loss(mv1, mv2)
    jensen_shannon_distance(mv1, mv2)
end

function _independent_loss(normals1, normals2)
    sum([jensen_shannon_distance(normals1[i], normals2[i]) for i in eachindex(normals1)])
end

_coord_vector(layers, c) = layers[c[1], c[2], :]

function _multivariate_full_dist(layers)
    numlayers =  size(layers,3)
    Σ = zeros(Float32, numlayers,numlayers)
    
    for i in 1:numlayers, j in i:numlayers
        if i != j 
            Σ[i,j] = cov(vec(layers[:,:,i]), vec(layers[:,:,j]))
            Σ[j,i] = Σ[i,j]
        else
            Σ[i,j] = std(vec(layers[:,:,i]))
        end
    end    
    μ = [mean(layers[:,:,i]) for i in 1:numlayers]
    MvNormal(μ,Σ)
end


function _multivariate_sampled_dist(sampled_layers)
    numlayers =  size(sampled_layers,2)
    Σ = zeros(Float32,numlayers,numlayers)
    
    for i in 1:numlayers, j in i:numlayers
        if i != j 
            Σ[i,j] = cov(sampled_layers[:,i],sampled_layers[:,j])
            Σ[j,i] = Σ[i,j]
        else
            Σ[i,j] = std(vec(sampled_layers[:,i]))
        end
    end    
    μ = [mean(sampled_layers[:,i]) for i in 1:numlayers]
    MvNormal(μ,Σ)
end

function _independent_full_dist(layers)
    dists = []
    for i in 1:size(layers,3)
        push!(dists, Normal(mean(layers[:, :, i]), std(layers[:,:,i])))
    end 
    dists
end

function _independent_sampled_dist(sampled_layers)
    numlayers = size(sampled_layers,2)
    [Normal(mean(sampled_layers[:,l]), std(sampled_layers[:,l])) for l in 1:numlayers]
end

function _generate!(
    coords::Vector{CartesianIndex},
    pool::Vector{CartesianIndex},
    sampler::KLSimulatedAnnealing,
    uncertainty,
) 
    layers = sampler.layers
    ndims(layers) <= 2 && throw(ArgumentError("KL-divergence Simulated Annealing needs more than one layer to work."))
    size(uncertainty) != (size(layers,1), size(layers,2)) && throw(DimensionMismatch("Layers are not the same dimension as uncertainty"))
    
    loss = sampler.multivariate ? _multivariate_loss : _independent_loss

    full_dist = sampler.multivariate ?  _multivariate_full_dist(layers) : _independent_full_dist(layers)
    
    np = sampler.numpoints
    coords .= pool 
    temps = temperature_schedule(sampler)
    
    log_freq = 50

    best_kl = zeros(Float32, (length(temps)*np+1))
    best_kl[1] = 1000.
    best_kl_overall = best_kl[begin]
 
    coord_sets = []

    sampled_layers = zeros(Float32, np, size(layers, 3))
    
    i = 2
    
    progbar = ProgressMeter.Progress(Int32(floor(length(best_kl)/log_freq)))
    accepts = 0
    for T in temps
        these_sampled_layers = deepcopy(sampled_layers)
        for p in 1:np
            proposed_coord = propose(sampler, coords[p])
            proposed_vec = _coord_vector(layers, proposed_coord)
            
            tmp_storage = zeros(Float32, size(layers,3))
            tmp_storage .= these_sampled_layers[p,:] 

            these_sampled_layers[p,:] .= proposed_vec
            sampled_dist = sampler.multivariate ? _multivariate_sampled_dist(these_sampled_layers) : _independent_sampled_dist(these_sampled_layers)

            KL_proposed = loss(full_dist, sampled_dist)
            if rand() <= move_probability(KL_proposed, best_kl[i - 1], T)
                accepts += 1
                if KL_proposed < best_kl_overall
                    best_kl[i] = KL_proposed
                    best_kl_overall = KL_proposed
                else
                    best_kl[i] = best_kl_overall
                end
            else
                these_sampled_layers[p,:]  .= tmp_storage 
                coords[p] = CartesianIndex(proposed_coord)
                best_kl[i] = best_kl[i - 1]
            end    

            if i % log_freq == 0
                ProgressMeter.next!(progbar; 
                showvalues = [
                        (Symbol("Best KL"), best_kl[i]),
                        (Symbol("Acceptance Ratio"), accepts/i)
                    ]
                )
                push!(coord_sets, copy(coords))
            end 


            i += 1
        end
    end
    return (coords, uncertainty), best_kl, coord_sets
end


layers = BiodiversityObservationNetworks.stack([rand(MidpointDisplacement(0.3), (1000,1000)) for _ in 1:2])
np= 100





klsa = KLSimulatedAnnealing(
    layers=layers, 
    T₀=10., 
    steps=10000,
    numpoints=np,
    cooling=GeometricCooling(0.999),
    multivariate=true
)

pool, unce = seed(BalancedAcceptance(numpoints=np), rand(1000,1000))
coords = similar(pool)

(_,_), kl, coord_set = _generate!(coords, pool, klsa, unce)

coord_set

klsa


coord_set[1]


a,b = MvNormal([0; 0], [1 0; 0 1]), MvNormal([2; 1],[4 0; 0 3])

jensen_shannon_distance(a,b)



using Plots

plot(eachindex(kl), log.(kl))
kl


# optim sandbo 


rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2




result = Optim.optimize(rosenbrock, zeros(2), SimulatedAnnealing(neighbor=a!), Optim.Options(iterations=10^6))

result = Optim.optimize(rosenbrock, zeros(2), NelderMead())

result = Optim.optimize(rosenbrock, zeros(2))


function a!(proposed, current)
    proposed .= (0.01rand() .+ current)
end


# end of optim sandbox






# -- foo
f()

#= MLP opt w/ zygote below, doesn't work 
sigmoid(x) = 1/(1+exp(-x))
softmax(x) = exp.(x) ./ sum(exp.(x))

mutable struct Linear
    W
    b
end
(l::Linear)(x) = sum(sigmoid.(l.W * x .+ l.b))


function get_sites!(model, layers, num_sites)
    output = zeros(Float32, size(layers)[1:2]...)
    I = CartesianIndices(Matrix(layers[:,:,1]))

    layer_buff = Zygote.Buffer(layers)
    x = Zygote.Buffer(output)
    for i in I
        v = model(layer_buff[i[1],i[2],:])
        x[i] = v
    end    
    y = copy(x)
    #CartesianIndices(y)[sortperm(vec(softmax(y)))[1:num_sites]]    
    #CartesianIndices(y)[sortperm(vec(softmax(y)))[1:num_sites]]
    sum(y)
end



function get_true_dists(layers)
    num_layers = size(layers)[3]
    [Normal(mean(layers[:,:,i]), var(layers[:,:,i])) for i in 1:num_layers]
end


# idea :
# is there any reason this can't consider MV gaussian across all layers,
# which would include covariance in KL div measure

function foo(layers, model, true_dists, num_epochs=20, num_sites=5)    
    sites = get_sites!(model, layers, num_sites)
    vals = [layers[s[1], s[2],1] for s in sites]
    kldivergence(true_dists[1], Normal(mean(vals), var(vals)))
    #sum([kldivergence(x...) for x in zip(true_dists, sample_dists)])
end

NUM_LAYERS = 5
model = Linear(rand(2,NUM_LAYERS), rand(2))
layers = BiodiversityObservationNetworks.stack([rand(MidpointDisplacement(0.5),100,200) for _ in 1:NUM_LAYERS])
true_dists = get_true_dists(layers)



get_sites!(model, layers,10)
Zygote.gradient(foo, layers, model, true_dists)


# brute force opt
emptmat = zeros(size(layers)[1:2]...)
nrounds = 1000
best_points, _ = seed(BalancedAcceptance(numpoints=30), emptmat)
best_kldiv = 1000.

vals = [layers[s[1], s[2],1] for s in best_points]

for i in 1:nrounds
    kldiv = 0.
    proposal, _ = seed(BalancedAcceptance(), emptmat)
    for l in 1:NUM_LAYERS
        vals = [layers[s[1], s[2],l] for s in proposal]
        k = kldivergence(true_dists[l], Normal(mean(vals), var(vals)))
        kldiv = kldiv + k
    end
    if kldiv < best_kldiv
        best_points = proposal
        best_kldiv = kldiv
    end
    @info best_kldiv
end 

using Plots

h1 = histogram(vec(layers[:,:,1]))
h2= histogram(vals)

plot(h1,h2)

scatter!([x[1] for x in best_points], [x[2] for x in best_points])

=#