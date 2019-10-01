# ### Load modules and include functions
include("../src/model_includes.jl")

#path names
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
results=dirname(pathof(PhytoAgentModel))*"/../results/"

# ### Read input files

#nTime = 1440 # number of time steps
nTime = 10 # number of time steps
ΔT = 3600 # time step: 3600 for 1 hour
temp,IR = read_input(samples*"T_IR.csv",trunc(Int,nTime*ΔT/3600));

RunParams=Dict("OutputResults"=>false,
               "Dims"=>3,
               "SaveTests"=>false,
               "NutOutputChoice"=>0)


g=load(samples*"grid.jld", "grid")
store_vel=load(samples*"uvw.jld", "uvw")

# ### Various model parameters

N = 100000   # Number of initial individuals of each species
Nsp = 2     # Number of species
Nn = Int(1e15)   # Number of cells one super-agent represents
B=setup_agents(N,Nsp,Cquota,Nn,1.0,0.25,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nut = [2.0, 0.05, 20.0, 0.0, 0.0, 0.0] #DIC, DIN, DOC, DON, POC, PON, mmol/m3
nutrients= setup_nutrients(g,nut)
if RunParams["NutOutputChoice"] == 2
    DIC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DIN = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DOC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    DON = zeros(g.Nx, g.Ny, g.Nz, nTime)
    POC = zeros(g.Nx, g.Ny, g.Nz, nTime)
    PON = zeros(g.Nx, g.Ny, g.Nz, nTime)
end
remin = remineralization(kDOC,kDON,kPOC,kPON);

# ### the main loop

for t in 1:nTime
    phyts_a = copy(B[t]) # read data from last time step
    phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
    velᵇ=store_vel[t]
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
    F = compute_nut_biochem(nutrients, remin)
    RunParams["Dims"] == 0 ? nothing : gtr = compute_source_term(nutrients, velᵇ, g, F)
    nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
    RunParams["OutputResults"] ? write_nut_cons(g, gtr, nutₜ, velᵇ, agent_num, t, death_ct, graz_ct, dvid_ct) : nothing
    if RunParams["NutOutputChoice"] == 1
        write_nut_nc_each_step(g, nutₜ, t) 
    elseif RunParams["NutOutputChoice"] == 2
        DIC[:,:,:,t] = nutₜ.DIC
        DIN[:,:,:,t] = nutₜ.DIN
        DOC[:,:,:,t] = nutₜ.DOC
        DON[:,:,:,t] = nutₜ.DON
        POC[:,:,:,t] = nutₜ.POC
        PON[:,:,:,t] = nutₜ.PON
    else
        nothing
    end
    global nutrients = nutₜ;
end

# ### post-processing steps
if RunParams["NutOutputChoice"] == 2
    write_nut_nc_alltime(g, DIC, DIN, DOC, DON, POC, PON, nTime)
end

B1 = []; B2 = [];
for i in 1:size(B,1)
    sort_species(B[i], B1, 1)
    sort_species(B[i], B2, 2)
end

if RunParams["Dims"] == 2 | RunParams["Dims"] == 3
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

if RunParams["Dims"] == 0 | RunParams["Dims"] == 2
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
RunParams["SaveTests"] ? CSV.write(results*"testB1B2_3D.csv",testB1B2(B1,B2)) : nothing
