# Colors for each type 
const _C_SAMPLER = Crayon(foreground = :blue, bold = true)
const _C_BON     = Crayon(foreground = :green, bold = true)
const _C_POOL    = Crayon(foreground = :magenta, bold = true)
const _C_UTIL    = Crayon(foreground = :yellow, bold = true)
const _C_DIM     = Crayon(foreground = :dark_gray)
const _C_RESET   = Crayon(reset = true) # used to turn of previous crayon to go back to default

# print text in a color
_c(io, crayon, text) = print(io, crayon, text, _C_RESET)

# Print sampler with fields in-line
function Base.show(io::IO, s::BONSampler)
    _c(io, _C_SAMPLER, nameof(typeof(s)))
    print(io, "(")
    fnames = fieldnames(typeof(s))
    for (i, f) in enumerate(fnames)
        i > 1 && print(io, ", ")
        _c(io, _C_DIM, "$f=")
        print(io, getfield(s, f))
    end
    print(io, ")")
end

# Print pool with number of sites and features
function Base.show(io::IO, cp::CandidatePool{K}) where {K}
    _c(io, _C_POOL, "CandidatePool")
    nf = ismissing(cp.features) ? 0 : size(cp.features, 1)
    print(io, "(", cp.n, " sites")
    nf > 0 && print(io, ", $(nf) features")
    print(io, ")")
end

# Print BON with number of sites and sampler used
function Base.show(io::IO, bon::BiodiversityObservationNetwork{K}) where {K}
    _c(io, _C_BON, "BiodiversityObservationNetwork")
    print(io, "(", length(bon.sites), " sites via ")
    _c(io, _C_SAMPLER, string(nameof(typeof(bon.sampler))))
    print(io, ")")
end

# Show metric with fields (if available)
function Base.show(io::IO, m::SamplingMetric)
    _c(io, _C_UTIL, nameof(typeof(m)))
    print(io, "(")
    fnames = fieldnames(typeof(m))
    for (i, f) in enumerate(fnames)
        i > 1 && print(io, ", ")
        _c(io, _C_DIM, "$f=")
        print(io, getfield(m, f))
    end
    print(io, ")")
end
