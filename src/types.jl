"""
    abstract type BONSeeder end

A `BONSeeder` is an algorithm for proposing sampling locations using a raster of
weights, represented as a matrix, in each cell.
"""
abstract type BONSeeder end

"""
    abstract type BONRefiner end 

A `BONRefiner` is an algorithm for proposing sampling locations by _refining_ a
set of candidate points to a smaller set of 'best' points.
"""
abstract type BONRefiner end

"""
    BONSampler

A union of the abstract types `BONSeeder` and `BONRefiner`. Both types return a
tuple with the coordinates as a vector of `CartesianIndex`, and the weight
matrix as a `Matrix` of `AbstractFloat`, in that order.
"""
const BONSampler = Union{BONSeeder, BONRefiner}


numsites(sampler::BONSampler) = sampler.numsites