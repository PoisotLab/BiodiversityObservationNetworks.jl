abstract type BONException <: Base.Exception end
abstract type SeederException <: BONException end

Base.showerror(io::IO, e::E) where {E <: BONException} =
    tprint(
        io,
        "{bold red}$(supertype(E)){/bold red} |  {bold yellow}$(e.message){/bold yellow}\n",
    )

function _check_arguments(sampler::S) where {S <: Union{BONSeeder, BONRefiner}}
    sampler.numsites > 1 || throw(TooFewSites())
    return nothing
end

@kwdef struct TooFewSites <: BONException
    message = "Number of sites to select must be at least two."
end
function check(::Type{TooFewSites}, sampler)
    sampler.numsites > 1 || throw(TooFewSites())
    return nothing
end

@kwdef struct TooManySites <: BONException
    message = "Cannot select more sites than there are candidates."
end
function check(::Type{TooManySites}, sampler)
    sampler.numsites <= maxsites(sampler) || throw(TooManySites())
    return nothing
end
