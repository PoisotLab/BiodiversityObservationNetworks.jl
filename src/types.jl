@kwdef mutable struct BiodiversityObservationNetwork{IT<:Integer,M<:AbstractMatrix}
    numobservatories::IT = 50
    coordinates::M = missing
end

abstract type BONSeeder end
abstract type BONRefiner end

BONSampler = Union{BONSeeder,BONRefiner}

function _checknumpoints(sampler::T) where {T<:BONSampler}
    getfield(sampler, :numpoints) > 0 || throw(ArgumentError("numpoints boooooo"))
    return nothing
end