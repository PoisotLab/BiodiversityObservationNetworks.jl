
function sample(alg::BONSampler) 
    _sample!(
        _allocate_sites(numsites(alg)), 
        _default_pool(alg),
        alg
    )
end


function sample(alg::BONSampler, l::L) where L<:Layer
    _sample!(
        _allocate_sites(numsites(alg)), 
        l,
        alg
    )
end


function sample(alg::BONSampler, candidates::C) where C<:Sites
    _sample!(
        _allocate_sites(numsites(alg)), 
        candidates,
        alg
    )
end
