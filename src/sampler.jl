"""
    BONSampler

Abstract supertype for all spatial sampling algorithms. Each concrete sampler
is a struct with at minimum an `n::Int` field specifying the desired
sample size.

Implement `_sample(rng, sampler, candidatepool)` to add a new algorithm.
"""
abstract type BONSampler end

# ========================================================================
# Sampler traits
# ========================================================================

"""Whether the sampler uses custom inclusion probabilities."""
supports_inclusion(::BONSampler) = false

"""Whether the sampler can support features (auxiliary variables associated with each site)."""
supports_features(::BONSampler) = false

"""Whether the sampler requires features (auxiliary variables associated with each site) to be present."""
requires_features(::BONSampler) = false

"""Whether the sampler guarantees exactly `n` selected sites."""
guarantees_exact_n(::BONSampler) = false