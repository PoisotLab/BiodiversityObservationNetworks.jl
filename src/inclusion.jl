"""
    get_uniform_inclusion(sampler, domain)

Construct a uniform inclusion surface for `domain` such that the sum of
inclusion probabilities equals `sampler.num_nodes`.

Returns a `RasterDomain` when given a raster-like domain.
"""
function get_uniform_inclusion(sampler, domain::T) where T<:Union{RasterDomain, RasterStack}
    pool = getpool(domain)
    inclusion = zeros(size(domain))
    for p in pool
        inclusion[p] = sampler.num_nodes / length(pool)
    end
    return RasterDomain(inclusion, domain.pool)
end

"""
    get_uniform_inclusion(sampler, bon::BiodiversityObservationNetwork)

Construct a uniform inclusion vector for a BON-like domain.
"""
function get_uniform_inclusion(sampler, domain::BiodiversityObservationNetwork)
    # TODO: this breaks if non-cartidx 
    pool = getpool(domain)
    inclusion = zeros(size(domain))
    for p in eachindex(pool)
        inclusion[p] = sampler.num_nodes / length(pool)
    end
    return inclusion
end
"""
    _redistribute_masked_weight!(domain, inclusion::RasterDomain{<:SDMLayer})

Redistribute any probability mass assigned to masked-out cells uniformly across
valid indices in `domain`. Mutates `inclusion` in-place.
"""
function _redistribute_masked_weight!(domain, inclusion::RasterDomain{<:SDMLayer})
    valid_idx = getpool(domain)

    weight_to_redistribute = sum(inclusion) - sum(inclusion.data.grid[valid_idx])
    weight_per_idx = weight_to_redistribute / length(valid_idx)

    inclusion.data.grid[findall(.!domain.pool)] .= 0
    inclusion.data.grid[valid_idx] .+= weight_per_idx
end

"""
    _redistribute_masked_weight!(domain, inclusion::RasterDomain{<:Matrix})

Matrix-backed variant of masked weight redistribution.
"""
function _redistribute_masked_weight!(domain, inclusion::RasterDomain{<:Matrix})
    valid_idx = getpool(domain)
    weight_to_redistribute = sum(inclusion) - sum(inclusion.data[valid_idx])
    weight_per_idx = weight_to_redistribute / length(valid_idx)

    inclusion.data[findall(.!domain.pool)] .= 0
    inclusion.data[valid_idx] .+= weight_per_idx
end

# ------------------------------------------
#  Inclusion normalization 
#
# ------------------------------------------
"""
    normalize_inclusion!(sampler, inclusion)

Scale `inclusion` so that its total mass equals `sampler.num_nodes`.
Overloads exist for SDMLayer-backed rasters, matrix-backed rasters, and vectors.
Mutates `inclusion` in-place.
"""
function normalize_inclusion!(sampler, inclusion::RasterDomain{<:SDMLayer})
    inclusion.data.grid .= sampler.num_nodes .* (inclusion.data.grid ./ sum(inclusion))
end 

function normalize_inclusion!(sampler, inclusion::RasterDomain{<:Matrix})
    inclusion.data .= sampler.num_nodes .* (inclusion.data ./ sum(inclusion))
end 

function normalize_inclusion!(sampler, inclusion::Vector)
    inclusion .= sampler.num_nodes .* (inclusion ./ sum(inclusion))
end 

# ------------------------------------------
#  Inclusion conversion  
#
# ------------------------------------------
"""
    convert_inclusion(sampler, domain, inclusion; kwargs...)

Coerce an `inclusion` specification to a representation compatible with `domain`.
Also ensures totals sum to `sampler.num_nodes` (renormalizing if needed) and
redistributes probability assigned to masked-out cells back to valid cells.
Returns `nothing` when `inclusion` is `nothing`.
"""
convert_inclusion(sampler, domain, ::Nothing; kwargs...) = nothing

function convert_inclusion(
    sampler, 
    domain::RasterDomain{<:SDMLayer}, 
    inclusion::AbstractMatrix; 
    kwargs...
) 
    x,y = extent(domain)
    rd = RasterDomain(SDMLayer(Matrix(inclusion); x=x, y=y, crs=domain.data.crs), domain.pool)
    convert_inclusion(sampler, domain, rd; kwargs...)
end 

function convert_inclusion(
    sampler, 
    domain::RasterDomain{<:Matrix}, 
    inclusion::AbstractMatrix; 
    kwargs...
) 
    rd = RasterDomain(Matrix(inclusion), domain.pool)
    convert_inclusion(sampler, domain, rd; kwargs...)
end 

function convert_inclusion(
    sampler, 
    domain::RasterDomain{<:SDMLayer},
    inclusion::SDMLayer; 
    kwargs...
)
    convert_inclusion(sampler, domain, RasterDomain(inclusion, domain.pool); kwargs...)
end 

function convert_inclusion(
    sampler, 
    domain::BiodiversityObservationNetwork,
    inclusion::Vector; 
    kwargs...
)
    if !isapprox(sum(inclusion), sampler.num_nodes)
        @warn "Inclusion probabilities do not sum to sampler's number of nodes. Interpreting them as weights"
        
        normalize_inclusion!(sampler, inclusion)
    end
    return inclusion 
    #convert_inclusion(sampler, domain, inclusion; kwargs...)
end 


# TODO: mask is unnecessary as an arg if the domain.pool is adjusted 


function convert_inclusion(
    sampler, 
    domain::RasterDomain, 
    inclusion::RasterDomain; 
    kwargs...
)
    size(domain) == size(inclusion) || throw(ArgumentError("Domain and Inclusion probabilities must be same size."))
    extent(domain) == extent(inclusion) || throw(ArgumentError("Domain and Inclusion probabilities must have same extent."))
    crs(domain) == crs(inclusion) || throw(ArgumentError("Domain and Inclusion probabilities must have CRS."))

    if !isapprox(sum(inclusion), sampler.num_nodes)
        @warn "Inclusion probabilities do not sum to sampler's number of nodes. Interpreting them as weights"
        
        normalize_inclusion!(sampler, inclusion)
    end

    _redistribute_masked_weight!(domain, inclusion)

    return inclusion
end