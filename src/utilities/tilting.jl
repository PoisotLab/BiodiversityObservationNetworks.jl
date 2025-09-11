"""
    tilt(layer, α)

Performs logistic-exponential tilting on the on a layer with scaling factor α. This is useful for adjusting inclusion probabilities to scale more toward (α > 1) or away from (α < 1) extreme values. 
"""
function tilt(layer, α)
    exp.(α * layer) ./ (1 .+ exp.(layer))
end