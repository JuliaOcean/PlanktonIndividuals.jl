# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# ### Load modules and include functions
include("model_includes.jl")

#path names
samples=dirname(pathof(PhytoAgentModel))*"/../samples/"
results=dirname(pathof(PhytoAgentModel))*"/../results/"

# ### Read input files

#nTime = 1440 # number of time steps
nTime = 10 # number of time steps
ΔT = 3600 # time step: 3600 for 1 hour
temp,IR = read_input(samples*"T_IR.csv",trunc(Int,nTime*ΔT/3600));

RunParams=Dict("OutputResults"=>false,"GridChoice"=>2,"VelChoice"=>2,
"SaveGrid"=>false,"SaveVel"=>false,"SaveTests"=>false)

function load_grid_0354(path)
# grid selected : 
  fieldroot = path*"/run.0354/";
  g = grid_offline(fieldroot);
end

RunParams["GridChoice"]==1 ? g=load_grid_0354(samples) : g=load(samples*"grid.jld", "grid")
RunParams["SaveGrid"] ? save(results*"grid.jld", "grid", g) : nothing

# ### remove old result files and create `results/` if needed
RunParams["OutputResults"] ? PrepRunDir(results) : nothing

# ### deal with time steps of offline velocity fields

itvalLo = 144;
itvalHi = 687888;
itList = collect(itvalLo:144:itvalHi);
tN = 3336; # starting time

vfroot = samples*"run.0354/offline-0604/" # directory of velocity fields
store_vel=[] #for storing and saving velocities when RunParams["SaveVel"]
RunParams["VelChoice"]==2 ? store_vel=load(samples*"uvw.jld", "uvw") : nothing

# ### Various model parameters

N = 100000   # Number of initial individuals of each species
Nsp = 2     # Number of species
Nn = Int(1e15)   # Number of cells one super-agent represents
B=setup_agents(N,Cquota,Nn,1.0,0.25,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nut = [2.0, 0.05, 20.0, 0.0, 0.0, 0.0] #DIC, DIN, DOC, DON, POC, PON, mmol/m3
nutrients= setup_nutrients(g,nut)
remin = remineralization(kDOC,kDON,kPOC,kPON);

# ### the main loop

for t in 1:nTime
    phyts_a = copy(B[t]) # read data from last time step
    phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
    if RunParams["VelChoice"]==1
        velᵇ = read_offline_vels(vfroot,itList,tN,trunc(Int,t*ΔT/3600))
    else
        velᵇ=store_vel[t]
    end
    global store_vel; RunParams["SaveVel"] ? store_vel=push!(store_vel,velᵇ) : nothing
    velᵈ = double_grid(velᵇ,g)
    agent_move(phyts_b,velᵈ,g,ΔT)
    push!(B,phyts_b)
    write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,output)
    agent_num = size(phyts_b,1)
    F = compute_nut_biochem(nutrients, remin)
    gtr = compute_source_term(nutrients, velᵇ, g, F)
    nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
    RunParams["OutputResults"] ? write_nut_nc(g, nutₜ, t) : nothing
    RunParams["OutputResults"] ? write_nut_cons(g, gtr, nutₜ, velᵇ, agent_num, t) : nothing
    global nutrients = nutₜ;
end

# ### post-processing steps

B1 = []; B2 = [];
for i in 1:size(B,1)
    sort_species(B[i], B1, B2)
end

HD1 = []; HD2 = [];
for i in 1:size(B,1)
    HD_1 = count_horizontal_num(B1[i],g);
    push!(HD1,HD_1)
    HD_2 = count_horizontal_num(B2[i],g);
    push!(HD2,HD_2)
end

for i in 1:size(B,1)
    convert_coordinates(B1[i],g) # convert grids to lon, lat and depth
    convert_coordinates(B2[i],g) # convert grids to lon, lat and depth
end

# ### more post-processing steps

VD1 = []; VD2 = [];
for i in 1:size(B,1)
    VD_1 = count_vertical_num(B1[i]);
    push!(VD1,VD_1)
    VD_2 = count_vertical_num(B2[i]);
    push!(VD2,VD_2)
end

output1, output2 = compute_mean_species(B1, B2, nTime);

# ### save to file

if RunParams["OutputResults"]
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
    RunParams["SaveVel"] ? save(results*"uvw.jld", "uvw", store_vel) : nothing
    RunParams["SaveTests"] ? CSV.write(results*"testB1B2.csv",testB1B2(B1,B2)) : nothing
end #if RunParams["OutputResults"]
