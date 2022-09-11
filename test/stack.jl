module BONTestStack
    using BiodiversityObservationNetworks
    using SimpleSDMLayers
    using Test

    nl = 10
    layers = [rand(50,50) for i in 1:nl]

    @test typeof(stack(layers)) <: Array{T,3} where T
    @test size(stack(layers),3) == nl

    bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7)
    temp, precip, seasonality, elevation = 
        convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 7; bbox...)),
        convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 12; bbox...)),
        convert(Float32, SimpleSDMPredictor(WorldClim, BioClim, 4; bbox...)),
        convert(Float32, SimpleSDMPredictor(WorldClim, Elevation; bbox...))
    
    layers = [temp,precip,seasonality,elevation] 
    @test typeof(stack(layers)) <: Array{T,3} where T
end 