abstract type BONException <: Base.Exception end
abstract type SeederException <: BONException end

Base.showerror(io::IO, e::E) where {E <: BONException} =
    tprint(
        io,
        "{bold red}$(supertype(E)){/bold red} |  {bold yellow}$(e.message){/bold yellow}\n",
    )

function _check_arguments(sampler::S) where {S <: Union{BONSeeder, BONRefiner}}
    return sampler.numsites > 1 || throw(TooFewSites(sampler.numsites))
end

@kwdef struct TooFewSites <: BONException
    message = "Number of sites to select must be at least two."
end
function check(TooFewSites, sampler)
    return sampler.numsites > 1 || throw(TooFewSites())
end

@kwdef struct TooManySites <: BONException
    message = "Cannot select more sites than there are candidates."
end
function check(TooManySites, sampler, max_sites)
    return sampler.numsites <= max_sites || throw(TooManySites())
end
