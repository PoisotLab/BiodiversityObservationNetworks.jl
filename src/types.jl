Base.@kwdef mutable struct BiodiversityObservationNetwork{IT<:Integer,M<:AbstractMatrix}
    numobservatories::IT = 50
    coordinates::M = missing
end

abstract type BONSeeder end
abstract type BONRefiner end

BONSampler = Union{BONSeeder,BONRefiner}
