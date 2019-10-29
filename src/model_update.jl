nTime = 10 # number of time steps
ΔT = 3600 # time step: 3600 for 1 hour
N = 100000   # Number of initial individuals of each species
Nsp = 2     # Number of species
Nn = Int(1e15)   # Number of cells one super-agent represents

samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
results=dirname(pathof(PhytoAgentModel))*"/../results/"

tmp = YAML.load_file(results*"params.yaml")
parameters = update_params(param_default, tmp)

#Define Grid Selection
Grid_sel = Dict("Nx"=>[551,560],"Ny"=>[1501,1510],"Nz"=>[1,40])

function load_grid_0354(path)
  fieldroot = path*"grid/run.0354/";
  g = grid_offline(fieldroot,Grid_sel);
end
GridOfflineOpt = Dict("gridpath" => dirname(pathof(PhytoAgentModel))*"/../samples/grid/run.0354/",
                      "GridSel" => Grid_sel)

RunOption.GridChoice ? g=load(samples*"grid.jld", "grid") : g=load_grid_0354(samples)
RunOption.SaveGrid ? save(results*"grid.jld", "grid", g) : nothing

# ### remove old result files and create `results/` if needed
RunOption.OutputResults ? PrepRunDir(results) : nothing

VelOfflineOpt = Dict("velpath" => dirname(pathof(PhytoAgentModel))*"/../samples/grid/run.0354/offline-0604",
                     "itList" => collect(144:144:687888), "tN" => 3336,
                     "GridSel" => Grid_sel)


vfroot = samples*"grid/run.0354/offline-0604" # directory of velocity fields
store_vel=[] #for storing and saving velocities when RunParams["SaveVel"]
RunOption.VelChoice ? store_vel=load(samples*"uvw.jld", "uvw") : nothing

# ### Various model parameters

B=setup_agents(N,Nsp,Cquota,Nn,1.0,0.25,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nut = [2.0, 0.05, 20.0, 0.0, 0.0, 0.0] #DIC, DIN, DOC, DON, POC, PON, mmol/m3
nutrients= setup_nutrients(g,nut)
gtr_0D = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
if RunOption.NutOutputChoice == false
    DIC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DIN = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DOC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DON = zeros(g.Nx, g.Ny, g.Nz, nTime)
    POC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    PON = zeros(g.Nx, g.Ny, g.Nz, nTime)
end

# ### the main loop

for t in 1:nTime
    phyts_a = copy(B[t]) # read data from last time step
    phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
    if RunOption.VelChoice == false
        velᵇ = read_offline_vels(vfroot,Grid_sel,itList,tN,trunc(Int,t*ΔT/3600))
    else
        velᵇ=store_vel[t]
    end
    global store_vel; RunOption.SaveVel ? store_vel=push!(store_vel,velᵇ) : nothing
    if g.Nx > 1 & g.Ny > 1 
        velᵈ = double_grid(velᵇ,g)
        agent_move(phyts_b,velᵈ,g,ΔT)
    elseif g.Nx == 1 & g.Ny == 1 & g.Nz > 1
        agent_move(phyts_b,velᵇ,g,ΔT) # for 1D only, use big grid velocities
    elseif g.Nx == 1 & g.Ny == 1 & g.Nz == 1
        nothing #for 0D only
    end
    push!(B,phyts_b)
    write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,output)
    agent_num = size(phyts_b,1)
    F = compute_nut_biochem(nutrients)
    RunParam.Dim == 0 ? gtr = gtr_0D : gtr = compute_source_term(nutrients, velᵇ, g, F)
    nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
    if RunOption.OutputResults
        write_nut_cons(g, gtr, nutₜ, velᵇ, agent_num, t, death_ct, graz_ct, dvid_ct)
        if RunOption.NutOutputChoice
            write_nut_nc_each_step(g, nutₜ, t) 
        else
            DIC[:,:,:,t] = nutₜ.DIC
            DIN[:,:,:,t] = nutₜ.DIN
            DOC[:,:,:,t] = nutₜ.DOC
            DON[:,:,:,t] = nutₜ.DON
            POC[:,:,:,t] = nutₜ.POC
            PON[:,:,:,t] = nutₜ.PON
        end
    else
        nothing
    end
    global nutrients = nutₜ;
end

# ### post-processing steps
if RunOption.OutputChoice == true & RunOption.NutOutputChoice == false
    write_nut_nc_alltime(g, DIC, DIN, DOC, DON, POC, PON, nTime)
end

B1 = []; B2 = [];
for i in 1:size(B,1)
    sort_species(B[i], B1, 1)
    sort_species(B[i], B2, 2)
end

if RunParam.Dim == 2 | RunParam.Dim == 3
    HD1 = []; HD2 = [];
    for i in 1:size(B,1)
        HD_1 = count_horizontal_num(B1[i],g);
        push!(HD1,HD_1)
        HD_2 = count_horizontal_num(B2[i],g);
        push!(HD2,HD_2)
    end
else
    nothing # for 1D or 0D
end

for i in 1:size(B,1)
    convert_coordinates(B1[i],g) # convert grids to lon, lat and depth
    convert_coordinates(B2[i],g) # convert grids to lon, lat and depth
end

# ### more post-processing steps
if RunParam.Dim == 0 | RunParam.Dim == 2
    nothing
else
    VD1 = []; VD2 = [];
    for i in 1:size(B,1)
        VD_1 = count_vertical_num(B1[i]);
        push!(VD1,VD_1)
        VD_2 = count_vertical_num(B2[i]);
        push!(VD2,VD_2)
    end
end

output1 = compute_mean_species(B1, nTime);
output2 = compute_mean_species(B2, nTime);

# ### save to file

if RunOption.OutputResults
    open(results*"B1.bin", "w") do io
        serialize(io, B1)
    end
    open(results*"B2.bin", "w") do io
        serialize(io, B2)
    end
    open(results*"grid.bin", "w") do io
        serialize(io, g)
    end
    open(results*"output.bin", "w") do io
        serialize(io, output)
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
    open(results*"IR.bin", "w") do io
        serialize(io, IR)
    end
end #if RunOption.OutputResults
RunOption.SaveVel ? save(results*"uvw.jld", "uvw", store_vel) : nothing
RunOption.SaveTests ? CSV.write(results*"testB1B2.csv",testB1B2(B1,B2)) : nothing
