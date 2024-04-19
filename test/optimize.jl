module BONTestOptimize
    using BiodiversityObservationNetworks
    using Test
    using Zygote
    using StatsBase

 #=   # empty method returns missing
    @test ismissing(optimize())

    layers = rand(50,50,5)
    loss = (layers, W, α) -> StatsBase.entropy(squish(layers, W, α))/prod(size(layers[:,:,1]))

    W,α,lossvals = optimize(layers, loss)
    @test typeof(W) <: Matrix
    @test typeof(α) <: Vector
    @test typeof(lossvals) <: Vector

    @test_throws ArgumentError optimize(layers, loss; learningrate = -0.1)
    @test_throws ArgumentError optimize(layers, loss; numsteps = 0)
    @test_throws ArgumentError optimize(zeros(50,50), loss; numsteps = 0)
    @test_throws ArgumentError optimize(layers, 5; numsteps = 0)=#

    

end 