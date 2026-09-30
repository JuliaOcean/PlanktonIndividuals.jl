using Test
using PlanktonIndividuals
using PlanktonKernels: CPU, RectilinearGrid
using PlanktonKernels.Units

@testset "PlanktonIndividuals" begin
    @testset "Unit tests" begin
        include("model_test.jl")
        include("parameter_test.jl")
        include("output_test.jl")
    end
    @testset "Example tests" begin
        include("test_0D_iron_energy_mode.jl")
        include("test_0D_protein_mode.jl")
        include("test_1D_macro_molecular_mode.jl")
        include("test_2D_quota_mode.jl")
        include("test_3D_carbon_mode.jl")
    end
end
