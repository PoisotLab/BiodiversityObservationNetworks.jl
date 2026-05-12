"""
    BiodiversityObservationNetwork{K}

The output of [`sample`](@ref). Holds selected sites together with the sampler
and domain that produced them, enabling downstream evaluation without passing
the domain separately.

# Fields
- `sites::Vector{K}`: selected keys 
- `coordinates::Matrix` — `d x n_selected` positions
- `features::Union{Matrix, Missing}` — `p x n_selected` or `Missing`
- `inclusion::Vector` — inclusion weights of the selected sites
- `sampler::BONSampler` — the sampler that produced this result
"""
struct BiodiversityObservationNetwork{K}
    sites::Vector{K}
    coordinates::Matrix
    features::Union{Matrix, Missing}
    inclusion::Vector
    sampler::BONSampler
end


# Base overloads
Base.length(bon::BiodiversityObservationNetwork) = length(bon.sites)
Base.size(bon::BiodiversityObservationNetwork) = (length(bon.sites),)
Base.iterate(bon::BiodiversityObservationNetwork, state = 1) = iterate(bon.sites, state)
Base.getindex(bon::BiodiversityObservationNetwork, i::Int) = bon.sites[i]
Base.firstindex(::BiodiversityObservationNetwork) = 1
Base.lastindex(bon::BiodiversityObservationNetwork) = length(bon.sites)

