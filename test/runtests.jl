
using Test
using PhytoAgentModel

@testset "PhytoAgentModel tests:" begin

#run("ln -s ../samples .")
#run("ln -s ../src .")

include(dirname(pathof(PhytoAgentModel))*"/model_update.jl")
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
tst=testB1B2(B1,B2,samples*"testB1B2.csv")
tst2=[tst[!,:B1],tst[!,:B1ref],tst[!,:B2],tst[!,:B2ref]]

@test isapprox(tst[end,:B1],tst[end,:B1ref]; atol=1e-2)
@test isapprox(tst[end,:B2],tst[end,:B2ref]; atol=1e-2)

end #@testset "PhytoAgentModel tests:" begin
