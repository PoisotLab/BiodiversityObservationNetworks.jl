"""
    _standardize

Standardizes the values of a matrix of predictors across the entire population
`Xfull`, and a set of predictors associated with the sampled sites, `Xsampled`
by scaling each dimension of the predictors to [0,1].

`Xsampled` is standardized based on the minimum and maximum values of each
predictor across the population, so both matrices a return on the same scale.

*Arguments*:
- `Xfull`: an `n` x `d` matrix, where `n` is the size of the population, and `d`
  is the number of predictors
- `Xsampled`: an `m` x `d` matrix, where `m` < `n` is the size of the sample

"""
function _standardize(Xfull, Xsampled)
    Xsampled_std = zeros(size(Xsampled))
    Xfull_std = zeros(size(Xfull))
    for i in axes(Xfull, 1)
        mi, mx = extrema(Xfull[i,:])
        Xfull_std[i,:] .= (Xfull[i,:] .- mi) ./ (mx - mi)
        Xsampled_std[i,:] = (Xsampled[i,:] .- mi) ./ (mx - mi)
    end 
    return Xfull_std, Xsampled_std
end

""" 
    jensenshannon

The [Jensen-Shannon
Divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence) is
a method for measuring the distance between two probability distibutions.

This method provides a comparison between the distribution of environmental
variables in a [`RasterStack`](@ref) `layers` to the values of those variables at the
sites within a [`BiodiversityObservationNetwork`](@ref) `bon`. 
"""
function jensenshannon(
    layers::RasterStack, 
    bon::BiodiversityObservationNetwork
)
    cart_idx, Xfull = features(layers)
    Xsampled = layers[bon]
    
    Xf, Xs = _standardize(Xfull, Xsampled)
    linear_idx = [findfirst(isequal(i), cart_idx) for i in BONs._get_cartesian_idx(layers, bon)]

    Xhat = zeros(size(Xf))
    Xhat[:, linear_idx] .= Xs

    M = 0.5(Xhat + Xf)

    return 0.5StatsBase.kldivergence(Xf, M) + 0.5StatsBase.kldivergence(Xf, M)
end 
