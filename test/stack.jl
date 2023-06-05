module BONTestStack
    using BiodiversityObservationNetworks
    using SpeciesDistributionToolkit
    using Test

    nl = 10
    layers = [rand(50,50) for i in 1:nl]

    @test typeof(BiodiversityObservationNetworks.stack(layers)) <: Array{T,3} where T
    @test size(BiodiversityObservationNetworks.stack(layers),3) == nl

    bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7)
    temp, precip, elevation = 
        SimpleSDMPredictor(rand(50,50); bbox...),
        SimpleSDMPredictor(rand(50,50); bbox...),
        SimpleSDMPredictor(rand(50,50); bbox...)
    
    layers = [temp,precip,elevation] 
    @test typeof(BiodiversityObservationNetworks.stack(layers)) <: Array{T,3} where T
end 
