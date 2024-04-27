abstract type BONException <: Base.Exception end
abstract type SeederException <: BONException end

Base.showerror(io::IO, e::E) where {E <: BONException} =
    tprint(
        io,
        "{bold red}$(supertype(E)){/bold red} |  {bold yellow}$E{/bold yellow}\n" *
        _error_message(e),
    )

function _check_arguments(sampler::S) where {S <: BONSeeder}
    return sampler.numpoints > 1 || throw(TooFewSites(sampler.numpoints))
end

struct TooFewSites{I} <: BONException where {I <: Integer}
    provided_sites::I
end
_error_message(tfs::TooFewSites) =
    "Number of sampling sites provided was $(tfs.provided_sites), but the number of sites must be {bold}greater than {cyan}1{/cyan}{/bold}.\n"

struct TooManySites{I} <: BONException where {I <: Integer}
    provided_sites::I
    maximum_sites::I
end
_error_message(tms::TooManySites) =
    "Number of sampling sites provided was $(tms.provided_sites), which {bold}exceeds the maximum possible{/bold} number of sites, $(tms.maximum_sites).\n"
