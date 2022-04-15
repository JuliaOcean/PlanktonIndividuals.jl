using Test
using PlanktonIndividuals

# @testset "PlanktonIndividuals 3D tests:" begin
#     include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_3D.jl")
#     @test isapprox(TC,TCt; atol=1e2)
# end #@testset "PlanktonIndividuals 3D tests:" begin

# @testset "PlanktonIndividuals 2D tests:" begin
#     include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_2D.jl")
#     @test isapprox(TP,TPt; atol=1e1)
# end #@testset "PlanktonIndividuals 2D tests:" begin

# @testset "PlanktonIndividuals 1D tests:" begin
#     include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_1D.jl")
#     @test isapprox(TP,TPt; atol=1e1)
# end #@testset "PlanktonIndividuals 1D tests:" begin

# @testset "PlanktonIndividuals 0D tests:" begin
#     include(dirname(pathof(PlanktonIndividuals))*"/../test/model_test_0D.jl")
#     @test isapprox(TP,TPt; atol=1e1)
# end #@testset "PlanktonIndividuals 0D tests:" begin

@testset "PlanktonIndividuals" begin
    @testset "Unit tests" begin
        include("grid_test.jl")
        include("field_test.jl")
    end
    @testset "Example tests" begin
        include("example_3D_test.jl")
        include("example_2D_test.jl")
        include("example_1D_test.jl")
        include("example_0D_test.jl")
    end
end
