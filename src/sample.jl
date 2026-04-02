


# ========================================================================
# Internal pipeline
# ========================================================================

# Checks if BONSampler is parameterized properly and can function with the given CandidatePool
function validate(sampler::BONSampler, cpool::CandidatePool)
    sampler.n > 0 || throw(ArgumentError("Sample size must be positive, got $(sampler.n)"))
    sampler.n <= cpool.n || throw(ArgumentError(
        "Sample size ($(sampler.n)) exceeds number of candidates ($(cpool.n))"))
    if requires_features(sampler) && ismissing(cpool.features)
        throw(ArgumentError(
            "$(typeof(sampler)) requires features but the CandidatePool has none."))
    end
end

# Fallback _sample, called when a sampler's internal method is not implemented, or the extension package is not loaded.
function _sample(::AbstractRNG, sampler::BONSampler, ::CandidatePool)
    error(
        "No sampling method is available for $(typeof(sampler)). " *
        "Using this sampler may require an extension package — see `?$(typeof(sampler))`.",
    )
end

# The internal call made by the public-facing `sample` method.
function _run_sample(rng::AbstractRNG, sampler::BONSampler, cpool::CandidatePool)
    validate(sampler, cpool)
    indices = _sample(rng, sampler, cpool)
    
    selected_keys = cpool.keys[indices]
    selected_coords = cpool.coordinates[:, indices]
    selected_feats = ismissing(cpool.features) ? Missing() : cpool.features[:, indices]
    selected_inclusion = cpool.inclusion[indices] ./ sum(cpool.inclusion[indices])
    BiodiversityObservationNetwork(selected_keys, selected_coords, selected_feats, selected_inclusion, sampler)
end


# ========================================================================
# Public sample API
# ========================================================================

"""
    sample([rng::AbstractRNG,] sampler::BONSampler, domain; mask=missing, inclusion=missing)

Select sites from `domain` using `sampler`. Pass an explicit `rng` for reproducibility.

`domain` can be any type accepted by [`CandidatePool`](@ref) — a `Matrix`,
a `Vector` of matrices, a [`CandidatePool`](@ref), or a [`BiodiversityObservationNetwork`](@ref)
(for multi-stage sampling).

`CandidatePool` can also accept types from [`SpeciesDistributionToolkit`](https://poisotlab.github.io/SpeciesDistributionToolkit.jl/),
including [`SDMLayer`]s from the [`SimpleSDMLayers`] subpackage, which represent rasters 
with geospatial metadata, and vector data types from the [`SimpleSDMPolygons`] subpackage.
The support for these types are included as extensions, meaning their functionality is only
loaded once [`SpeciesDistributionToolkit`] (or one of the corresponding subpackages) is loaded.
For more information on this functionality, see [this]() TODO how-to. 
"""
function sample end

# Without RNG
function sample(sampler::BONSampler, domain; mask=missing, inclusion=missing)
    sample(Random.default_rng(), sampler, domain; mask, inclusion)
end

function sample(sampler::BONSampler, cpool::CandidatePool)
    sample(Random.default_rng(), sampler, cpool)
end

function sample(sampler::BONSampler, bon::BiodiversityObservationNetwork)
    sample(Random.default_rng(), sampler, bon)
end

# With explicit RNG
function sample(rng::AbstractRNG, sampler::BONSampler, domain; mask=missing, inclusion=missing)
    cpool = CandidatePool(domain; mask, inclusion)
    _run_sample(rng, sampler, cpool)
end

function sample(rng::AbstractRNG, sampler::BONSampler, cpool::CandidatePool)
    _run_sample(rng, sampler, cpool)
end

function sample(rng::AbstractRNG, sampler::BONSampler, result::BiodiversityObservationNetwork)
    _run_sample(rng, sampler, CandidatePool(result))
end




#= 
function check_args(domain, mask, inclusion)
    if domain isa Matrix && mask isa Union{PolygonDomain,SimpleSDMPolygons.AbstractGeometry}
        throw(ArgumentError("Cannot use a polygon to mask a Matrix domain"))
    end
end

function preprocess(sampler, domain, mask, inclusion)
    check_args(domain, mask, inclusion)

    domain = to_domain(domain)
    mask = convert_mask(domain, mask)
    mask!(domain, mask)

    inclusion = convert_inclusion(sampler, domain, inclusion)

    return domain, mask, inclusion
end

function postprocess(domain, selected, auxiliary)
    BiodiversityObservationNetwork(selected, auxiliary)
end

function postprocess(domain::RasterDomain{<:SDMLayer}, selected, auxiliary)
    Es, Ns = eastings(domain.data), northings(domain.data)
    BiodiversityObservationNetwork([(Es[s[2]], Ns[s[1]]) for s in selected], auxiliary)
end

function postprocess(domain::BiodiversityObservationNetwork, selected, auxiliary)
    BiodiversityObservationNetwork(domain.nodes[selected], auxiliary)
end


# ------------------------------------------
#  Sample (user entrypoint)
#
# ------------------------------------------

"""
    sample
""" 
function sample() end

function sample(sampler::BONSampler) 
    sample(sampler, SDMLayer(zeros(180,90)))
end


function sample(
    sampler::BONSampler,
    domain;
    mask = nothing,
    inclusion = nothing,
    kwargs...
)
    domain, mask, inclusion = preprocess(sampler, domain, mask, inclusion)
    result, auxiliary = _sample(sampler, domain; inclusion=inclusion, kwargs...)
    #@info "Sampler: $sampler"
    #@info "Domain: $domain"
    #@info "Result: $result"
    #@info "Aux: $auxiliary"
    #@info "\n\n"
    postprocess(domain, result, auxiliary)
end


# ------------------------------------------
#  Tests
#
# ------------------------------------------

@testitem "We can sample from a matrix" setup=[TestModule] begin
    mat = zeros(20, 30)
    bon = sample(SimpleRandom(), mat)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    bon = sample(SimpleRandom(), layer)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:Matrix}" setup=[TestModule] begin
    rd = RasterDomain(zeros(20, 30))
    bon = sample(SimpleRandom(), rd)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:SDMLayer}" setup=[TestModule] begin
    rd = RasterDomain(SDT.SDMLayer(zeros(20, 30)))
    bon = sample(SimpleRandom(), rd)
    @test bon isa BiodiversityObservationNetwork
end

# ------------------------------------------
#  Masking 
#
# ------------------------------------------
@testitem "We can sample from a Matrix with a Matrix Mask" setup=[TestModule] begin
    mat = zeros(20, 30)
    mask = rand(size(mat)...) .> 0.1
    bon = sample(SimpleRandom(), mat, mask=mask)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer with a Matrix Mask of the same size" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    mask = rand(size(layer)...) .> 0.1
    bon = sample(SimpleRandom(), layer, mask=mask)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer with an SDMLayer Mask" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    mask_layer = SDT.SDMLayer(zeros(20, 30))
    mask_layer.indices .= rand(size(layer)...) .> 0.1
    bon = sample(SimpleRandom(), layer, mask=mask_layer)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:Matrix} with a Matrix Mask" setup=[TestModule] begin
end



@testitem "We can sample from a RasterDomain{<:SDMLayer} with a Matrix Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a SDMLayer Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a RasterDomain{<:Matrix} Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a RasterDomain{<:SDMLayer} Mask" setup=[TestModule] begin

end
=#

