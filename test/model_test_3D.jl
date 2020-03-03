using PhytoAgentModel, JLD
#path names
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
#                   Dim output, NutOutput, GridChoice, Gridoff, VelChoice, Veloff, SaveGrid, SaveVel, Test
RunOption=RunOptions(3, false,  true,      true,       Dict(),  true,      Dict(), false,    false,   false)

g=load(samples*"grid.jld", "grid");

model = PA_Model(g, RunParam;
              nutrients = setup_nutrients(g,[2.0, 0.05, 20.0, 0.0, 0.0, 0.0])); #DIC, DIN, DOC, DON, POC, PON, mmol/m3

PA_ModelRun(model, RunParam, RunOption)

# ### post-processing steps
B1 = []; B2 = [];
for i in 1:size(model.individuals,1)
    sort_species(model.individuals[i], B1, 1)
    sort_species(model.individuals[i], B2, 2)
end

output1 = compute_mean_species(B1, RunParam.nTime);
output2 = compute_mean_species(B2, RunParam.nTime);

RunOption.SaveTest ? CSV.write(samples*"testB1B2_3D.csv",testB1B2(B1,B2)) : nothing
