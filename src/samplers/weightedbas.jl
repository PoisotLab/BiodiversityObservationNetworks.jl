"""
    WeightedBalancedAcceptance 

Weighted Balaced Acceptance Sampling 


Any spatial design could in principle be weighted w/ inclusion probability
adjusted by a raster. (For some sampers this would be reall dumb, e.g. Grid. But
for spatially balanced, it makes sense )
"""
Base.@kwdef struct WeightedBalancedAcceptance{I<:Integer, F<:Real} <: BONSampler
    number_of_nodes::I = 100
    grid_size::Tuple{I,I} = (250, 250)
    inclusion_scaling::F = 1.
end 
WeightedBalancedAcceptance(n::Integer; grid_size=(250, 250)) = WeightedBalancedAcceptance(n, grid_size)

__valid_geometries(::WeightedBalancedAcceptance) = (Raster)

function _zscore(raster::SDMLayer) 
    z_raster = copy(raster)
    valid_idx = raster.indices
    Z = StatsBase.fit(ZScoreTransform, z_raster.grid[valid_idx])

    z_raster.grid[valid_idx] .= StatsBase.transform(Z, z_raster.grid[valid_idx])
    return z_raster
end 

function _sample(sampler::WeightedBalancedAcceptance, raster::SDMLayer)
    α = sampler.inclusion_scaling
    N = sampler.number_of_nodes
    inclusion_probability(x) = exp(α*x)/(1+exp(α*x))

    z_raster = _zscore(raster)

    Es, Ns = eastings(raster), northings(raster)
    x_dim, y_dim = length(Es), length(Ns)

    # TODO: this is heavily readundant w/ normal BAS, should reduce code duplication
    seed = rand(Int.(1e0:1e7), 2)
    selected_points = Node[]
    ct = 0
    candct = 0
    while ct < N
        i, j = haltonvalue(seed[1] + candct, 2), haltonvalue(seed[2] + candct, 3)
        candct += 1
        candx, candy = convert.(Int, [ceil(y_dim * i), ceil(x_dim * j)])
        candidate = CartesianIndex(candx,candy)
        cand_value = z_raster[candidate]

        if !isnothing(cand_value) && rand() < inclusion_probability(cand_value)
            push!(selected_points, Node((Es[candidate[2]], Ns[candidate[1]])))
            ct += 1
        end
    end
    return BiodiversityObservationNetwork(selected_points)
end

# ---------------------------------------------------------------
# Tests
# ---------------------------------------------------------------

@testitem "We can use the default Weighted BAS constructor on a raster" begin
    raster = BiodiversityObservationNetworks.SpeciesDistributionToolkit.SDMLayer(rand(50,30))

    wbas = WeightedBalancedAcceptance()
    bon = sample(wbas, raster)
    
    @test bon isa BiodiversityObservationNetwork
    @test size(bon) == wbas.number_of_nodes
end


