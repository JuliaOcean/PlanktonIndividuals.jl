using Test
using PhytoAgentModel

@testset "PhytoAgentModel 3D tests:" begin
    include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_3D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PhytoAgentModel 3D tests:" begin

@testset "PhytoAgentModel 2D tests:" begin
    include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_2D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PhytoAgentModel 2D tests:" begin

@testset "PhytoAgentModel 1D tests:" begin
    include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_1D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PhytoAgentModel 1D tests:" begin

@testset "PhytoAgentModel 0D tests:" begin
    include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_0D.jl")
    @test isapprox(TP,TPt; atol=1e1)
end #@testset "PhytoAgentModel 0D tests:" begin
