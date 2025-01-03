module BONsMakieExt

@info "Loading Makie Extension for BiodiversityObservationNetworks.jl..."

@static if isdefined(Base, :get_extension)
    using Makie, BiodiversityObservationNetworks
else    
    using ..Makie, ..BiodiversityObservationNetworks
end

include(joinpath(@__DIR__, "recipes.jl"))

end