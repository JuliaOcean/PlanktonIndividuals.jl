using Test
using PlanktonIndividuals

@testset "PlanktonIndividuals 3D tests:" begin
    include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_3D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PlanktonIndividuals 3D tests:" begin

@testset "PlanktonIndividuals 2D tests:" begin
    include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_2D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PlanktonIndividuals 2D tests:" begin

@testset "PlanktonIndividuals 1D tests:" begin
    include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_1D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PlanktonIndividuals 1D tests:" begin

@testset "PlanktonIndividuals 0D tests:" begin
    include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_0D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PlanktonIndividuals 0D tests:" begin
