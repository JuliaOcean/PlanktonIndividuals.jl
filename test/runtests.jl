
using Test
using PhytoAgentModel

@testset "PhytoAgentModel 3D tests:" begin

#run("ln -s ../samples .")
#run("ln -s ../src .")

include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_3D.jl")
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
tst_3D=testB1B2(B1,B2,samples*"testB1B2_3D.csv")
tst2=[tst_3D[!,:B1],tst_3D[!,:B1ref],tst_3D[!,:B2],tst_3D[!,:B2ref]]

@test isapprox(tst_3D[end,:B1],tst_3D[end,:B1ref]; atol=1e-2)
@test isapprox(tst_3D[end,:B2],tst_3D[end,:B2ref]; atol=1e-2)

end #@testset "PhytoAgentModel 3D tests:" begin

# @testset "PhytoAgentModel 2D tests:" begin

# include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_2D.jl")
# samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
# tst_2D=testB1B2(B1,B2,samples*"testB1B2_2D.csv")
# tst2=[tst_2D[!,:B1],tst_2D[!,:B1ref],tst_2D[!,:B2],tst_2D[!,:B2ref]]

# @test isapprox(tst_2D[end,:B1],tst_2D[end,:B1ref]; atol=1e-2)
# @test isapprox(tst_2D[end,:B2],tst_2D[end,:B2ref]; atol=1e-2)

# end #@testset "PhytoAgentModel 2D tests:" begin

# @testset "PhytoAgentModel 1D tests:" begin

# include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_1D.jl")
# samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
# tst_1D=testB1B2(B1,B2,samples*"testB1B2_1D.csv")
# tst2=[tst_1D[!,:B1],tst_1D[!,:B1ref],tst_1D[!,:B2],tst_1D[!,:B2ref]]

# @test isapprox(tst_1D[end,:B1],tst_1D[end,:B1ref]; atol=1e-2)
# @test isapprox(tst_1D[end,:B2],tst_1D[end,:B2ref]; atol=1e-2)

# end #@testset "PhytoAgentModel 1D tests:" begin

# @testset "PhytoAgentModel 0D tests:" begin

# include(dirname(pathof(PhytoAgentModel))*"/../test/model_test_0D.jl")
# samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
# tst_0D=testB1B2(B1,B2,samples*"testB1B2_0D.csv")
# tst2=[tst_0D[!,:B1],tst_0D[!,:B1ref],tst_0D[!,:B2],tst_0D[!,:B2ref]]

# @test isapprox(tst_0D[end,:B1],tst_0D[end,:B1ref]; atol=1e-2)
# @test isapprox(tst_0D[end,:B2],tst_0D[end,:B2ref]; atol=1e-2)

# end #@testset "PhytoAgentModel 0D tests:" begin
