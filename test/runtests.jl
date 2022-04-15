using Test
using PlanktonIndividuals

@testset "PlanktonIndividuals" begin
    @testset "Unit tests" begin
        include("grid_test.jl")
        include("field_test.jl")
        include("parameter_test.jl")
        include("output_test.jl")
    end
    @testset "Example tests" begin
        include("example_3D_test.jl")
        include("example_2D_test.jl")
        include("example_1D_test.jl")
        include("example_0D_test.jl")
    end
end
