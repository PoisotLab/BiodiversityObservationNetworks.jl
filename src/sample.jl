function check_args(domain, mask, inclusion)
    if domain isa Matrix && mask isa Union{PolygonDomain,SimpleSDMPolygons.AbstractGeometry}
        throw(ArgumentError("Cannot use a polygon to mask a Matrix domain"))
    end
end

function preprocess(sampler, domain, mask, inclusion)
    check_args(domain, mask, inclusion)

    domain = to_domain(domain)
    mask = convert_mask(domain, mask)
    mask!(domain, mask)

    inclusion = convert_inclusion(sampler, domain, inclusion)

    return domain, mask, inclusion
end

function postprocess(domain, selected, auxiliary)
    BiodiversityObservationNetwork(selected, auxiliary)
end

function postprocess(domain::BiodiversityObservationNetwork, selected, auxiliary)
    BiodiversityObservationNetwork(domain.nodes[selected], auxiliary)
end


# ------------------------------------------
#  Sample (user entrypoint)
#
# ------------------------------------------

"""
    sample
""" 
function sample() end

function sample(sampler::BONSampler) 
    sample(sampler, SDMLayer(zeros(180,90)))
end


function sample(
    sampler::BONSampler,
    domain;
    mask = nothing,
    inclusion = nothing,
    kwargs...
)
    domain, mask, inclusion = preprocess(sampler, domain, mask, inclusion)
    result, auxiliary = _sample(sampler, domain; inclusion=inclusion, kwargs...)
    #@info "Sampler: $sampler"
    #@info "Domain: $domain"
    #@info "Result: $result"
    #@info "Aux: $auxiliary"
    #@info "\n\n"
    postprocess(domain, result, auxiliary)
end


# ------------------------------------------
#  Tests
#
# ------------------------------------------

@testitem "We can sample from a matrix" setup=[TestModule] begin
    mat = zeros(20, 30)
    bon = sample(SimpleRandom(), mat)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    bon = sample(SimpleRandom(), layer)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:Matrix}" setup=[TestModule] begin
    rd = RasterDomain(zeros(20, 30))
    bon = sample(SimpleRandom(), rd)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:SDMLayer}" setup=[TestModule] begin
    rd = RasterDomain(SDT.SDMLayer(zeros(20, 30)))
    bon = sample(SimpleRandom(), rd)
    @test bon isa BiodiversityObservationNetwork
end

# ------------------------------------------
#  Masking 
#
# ------------------------------------------
@testitem "We can sample from a Matrix with a Matrix Mask" setup=[TestModule] begin
    mat = zeros(20, 30)
    mask = rand(size(mat)...) .> 0.1
    bon = sample(SimpleRandom(), mat, mask=mask)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer with a Matrix Mask of the same size" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    mask = rand(size(layer)...) .> 0.1
    bon = sample(SimpleRandom(), layer, mask=mask)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a SDMLayer with an SDMLayer Mask" setup=[TestModule] begin
    layer = SDT.SDMLayer(zeros(20, 30))
    mask_layer = SDT.SDMLayer(zeros(20, 30))
    mask_layer.indices .= rand(size(layer)...) .> 0.1
    bon = sample(SimpleRandom(), layer, mask=mask_layer)
    @test bon isa BiodiversityObservationNetwork
end

@testitem "We can sample from a RasterDomain{<:Matrix} with a Matrix Mask" setup=[TestModule] begin
end



@testitem "We can sample from a RasterDomain{<:SDMLayer} with a Matrix Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a SDMLayer Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a RasterDomain{<:Matrix} Mask" setup=[TestModule] begin

end

@testitem "We can sample from a RasterDomain{<:SDMLayer} with a RasterDomain{<:SDMLayer} Mask" setup=[TestModule] begin

end


# ------------------------------------------
#  Custom Inclusion Probabilities 
#
# ------------------------------------------


# ------------------------------------------
#  Mask and Custom Inclusion Probabilities 
#
# ------------------------------------------



