using PhytoAgentModel, JLD, Serialization, YAML
#path names
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
results=dirname(pathof(PhytoAgentModel))*"/../results/"
#                   Dim output, NutOutput, GridChoice, Gridoff, VelChoice, Veloff, SaveGrid, SaveVel, Test
RunOption=RunOptions(3, true,   true,      false,      Dict(),  false,     Dict(), false,    false,   false)

tmp = YAML.load_file(results*"params.yaml")
parameters = update_params(param_default, tmp)

#Define Grid Selection
Grid_sel = Dict("Nx"=>[551,560],"Ny"=>[1501,1510],"Nz"=>[1,40])

RunOption.GridOfflineOpt = Dict("gridpath" => dirname(pathof(PhytoAgentModel))*"/../samples/grid/run.0354/",
                      "GridSel" => Grid_sel)



RunOption.GridChoice ? g=load(samples*"grid.jld", "grid") : g = grid_offline(RunOption.GridOfflineOpt);
RunOption.SaveGrid ? save(results*"grid.jld", "grid", g) : nothing

# ### remove old result files and create `results/` if needed
RunOption.OutputChoice ? PrepRunDir(results) : nothing

RunOption.VelOfflineOpt = Dict("velpath" => dirname(pathof(PhytoAgentModel))*"/../samples/grid/run.0354/offline-0604",
                     "itList" => collect(144:144:687888), "tN" => 3336,
                     "GridSel" => Grid_sel)

model = Model(g, RunParam;
              nutrients = setup_nutrients(g,[2.0, 0.05, 20.0, 0.0, 0.0, 0.0])) #DIC, DIN, DOC, DON, POC, PON, mmol/m3

ModelRun(model, RunParam, RunOption)

# ### post-processing steps
B1 = []; B2 = [];
for i in 1:size(model.individuals,1)
    sort_species(model.individuals[i], B1, 1)
    sort_species(model.individuals[i], B2, 2)
end

HD1 = []; HD2 = [];
for i in 1:size(model.individuals,1)
    HD_1 = count_horizontal_num(B1[i],g);
    push!(HD1,HD_1)
    HD_2 = count_horizontal_num(B2[i],g);
    push!(HD2,HD_2)
end

for i in 1:size(model.individuals,1)
    convert_coordinates(B1[i],g) # convert grids to lon, lat and depth
    convert_coordinates(B2[i],g) # convert grids to lon, lat and depth
end

output1 = compute_mean_species(B1, RunParam.nTime);
output2 = compute_mean_species(B2, RunParam.nTime);

VD1 = []; VD2 = [];
for i in 1:size(model.individuals,1)
    VD_1 = count_vertical_num(B1[i]);
    push!(VD1,VD_1)
    VD_2 = count_vertical_num(B2[i]);
    push!(VD2,VD_2)
end

RunOption.SaveTest ? CSV.write(results*"testB1B2_3D.csv",testB1B2(B1,B2)) : nothing

# ### save to file

if RunOption.OutputResults
    open(results*"B1.bin", "w") do io
        serialize(io, B1)
    end
    open(results*"B2.bin", "w") do io
        serialize(io, B2)
    end
    open(results*"grid.bin", "w") do io
        serialize(io, model.grid)
    end
    open(results*"output.bin", "w") do io
        serialize(io, model.output)
    end
    open(results*"output1.bin", "w") do io
        serialize(io, output1)
    end
    open(results*"output2.bin", "w") do io
        serialize(io, output2)
    end
    open(results*"VD1.bin", "w") do io
        serialize(io, VD1)
    end
    open(results*"VD2.bin", "w") do io
        serialize(io, VD2)
    end
    open(results*"HD1.bin", "w") do io
        serialize(io, HD1)
    end
    open(results*"HD2.bin", "w") do io
        serialize(io, HD2)
    end
end #if RunOption.OutputResults
RunOption.SaveVel ? save(results*"uvw.jld", "uvw", store_vel) : nothing
RunOption.SaveTest ? CSV.write(results*"testB1B2.csv",testB1B2(B1,B2)) : nothing
