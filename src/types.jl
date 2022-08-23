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

A `BONRefiner` is an algorithm for proposing sampling locations
using a set of points with values as the 
"""
abstract type BONRefiner end

"""
    BONSampler

A combination of the abstract types `BONSeeder` and `BONRefiner`.
"""
const BONSampler = Union{BONSeeder,BONRefiner}