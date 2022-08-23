Base.@kwdef mutable struct BiodiversityObservationNetwork{IT<:Integer,M<:AbstractMatrix}
    numobservatories::IT = 50
    coordinates::M = missing
end

"""
    abstract type BONSeeder end

A `BONSeeder` is an algorithm for proposing sampling locations
using a raster of weights in each cell.
"""
abstract type BONSeeder end

"""
    abstract type BONRefiner end 
"""
abstract type BONRefiner end

"""
    
"""
const BONSampler = Union{BONSeeder,BONRefiner}