# Need to write tests for:
# - squish
# - optimize
# - stack
# - integrations

module BONTestSquish
    using BiodiversityObservationNetworks
    using NeutralLandscapes
    using Distributions
    using Test

    function makelayers(nl) 
        nl = 5
        layers = zeros(50,50,nl)
        for l in 1:nl
            layers[:,:,l] .= rand(MidpointDisplacement(), 50, 50)
        end
        layers
    end 

    layers = makelayers(5)
    ntargs = 3
    α = rand(Dirichlet([1 for t in 1:ntargs]))
    W = rand(5, ntargs) 
    for i in 1:size(W,2)
        W[:,i] .= W[:,i] ./ sum(W[:,i])
    end
    @test typeof(squish(layers, W, α)) <: Matrix

    layers = makelayers(5)
    ntargs = 3
    α = rand(Dirichlet([1 for t in 1:ntargs]))
    W = rand(5, ntargs)
    @test_throws ArgumentError squish(layers, W, α)

    W = rand(4, ntargs)
    @test_throws ArgumentError squish(layers, W, α)

    
    α = [1., 2., 3.]
    @test_throws ArgumentError squish(layers, W, α)

    α = rand(Dirichlet([1 for t in 1:ntargs+1]))
    @test_throws ArgumentError squish(layers, W, α)

    layers = makelayers(5)
    W = rand(4, ntargs)
    α = rand(Dirichlet([1 for t in 1:ntargs]))
    @test_throws ArgumentError squish(layers, W, α)

end