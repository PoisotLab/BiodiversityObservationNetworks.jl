Base.@kwdef struct Grid{F<:Real} <: BONSampler
    longitude_spacing::F = 1. # in wgs84 coordinates
    latitude_spacing::F  = 1.
end 


function _generate_grid(sampler::Grid, domain) 
    x, y = GeoInterface.extent(domain)
    x_step, y_step = sampler.longitude_spacing, sampler.latitude_spacing
    BiodiversityObservationNetwork([Node((i,j)) for i in x[1]:x_step:x[2], j in y[1]:y_step:y[2] if GeometryOps.contains(domain, (i,j))])
end 


function _sample(sampler::Grid, domain::T) where T 
    if GeoInterface.isgeometry(domain)
        return _generate_grid(sampler, domain)
    elseif GeoInterface.israster(domain)
        return 
    end 
    @warn "Can't use Grid on a $T"
end 
