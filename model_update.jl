# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

using DataFrames, CSV, NetCDF, ProgressMeter
using Random
using Distributions
cd("/nobackup1b/users/zhenwu/ABPM_3D/")
include("parameters.jl")
include("model_setup.jl")
include("model_struct.jl")
include("phyt_process.jl")
include("utils.jl")
include("agent_div.jl")
include("dst3fl.jl")
include("nutrient_processes.jl")
include("flux_div_diffusion_operators.jl")
# Read input files
nTime = 600 # number of time steps
ΔT = 3600 # time step: 1h
temp,IR = read_input("T_IR.csv",nTime);

# grid selected : [41:1040,1201:2600]: 26.46N to 47.71N, 172.188W to 151.375W
fieldroot = "/nobackup1b/users/jahn/hinpac/grazsame3/run/run.0354/";
g = grid_offline(fieldroot);

# deal with time steps of offline velocityfields
itvalLo = 144;
itvalHi = 687888;
itList = collect(itvalLo:144:itvalHi);
tN = 4056; # starting time
vfroot = "/nobackup1b/users/jahn/hinpac/grazsame3/run/run.0354/offline-0604/"; # directory of velocity fields

N = 50000   # Number of initial individuals of each species
Nsp = 2     # Number of species
B=setup_agents(N,Cquota,1.1,0.18,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nut = [2.0, 1.0, 20.0, 2.0, 1.0, 1.0] #DIC, DIN, DOC, DON, POC, PON, mmol/m3
nut₀= setup_nutrients(g,nut)
#nut₀.DIN[100,100,1] = 1.0
CN = [nut₀]
remin = rem(kDOC,kDON,kPOC,kPON)
#nutrients = DataFrame(DIN=1.0e-5, DOC=0.0, DON=3.0e-5, POC=0.0, PON=0.0);# mmol

@showprogress 1 "Updating..." for t in 1:600
    phyts_a = copy(B[t]) # read data from last time step
    nutrients = CN[t]
    velᵇ = read_offline_vels(vfroot,itList,tN,t);
    velᵈ = double_grid(velᵇ,g)
    agent_move(phyts_a,velᵈ,g,ΔT) 
    phyts_b,dvid_ct,graz_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
    push!(B,phyts_b)
    write_output(t,phyts_b,dvid_ct,graz_ct,output)
    convert_coordinates(B[t],g) # convert grids to lon, lat and depth
    F = compute_nut_biochem(nutrients, remin)
    gtr = compute_source_term(nutrients, velᵇ, g, F)
    nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
    #println("t_DIN:"*string(sum(nutₜ.DIN .* g.V)+sum(nutₜ.DON .* g.V)+sum(nutₜ.PON .* g.V))) #
    push!(CN, nutₜ)
end

B1 = []; B2 = [];
@showprogress 1 "Computing..." for i in 1:600
    sort_species(B[i], B1, B2)
end

output1, output2 = compute_mean_species(B1, B2, nTime);

###########################################################
# to save all the agents of all time steps run code below #
###########################################################
using JLD, FileIO
f = JLD.jldopen("results/output.jld", "w", compress=true)
JLD.@write f B1
JLD.@write f B2
JLD.@write f output
JLD.@write f output1
JLD.@write f output2
JLD.@write f CN
JLD.@write f g
JLD.@write f IR
close(f)

using Plots, LaTeXStrings
gr()
p1 = plot(output1.Population, xlims=(0,720), ylims=(0.0,135000.0), title="Population", xticks = 0:24:720, rotation=45, label="Population sp1", legend=:topleft);
plot!(p1,output2.Population, xlims=(0,720), xticks = 0:24:720, rotation=45,label="Population sp2");
plot!(p1,output.dvid, xlims=(0,720), xticks = 0:24:720, rotation=45,label="Cell Divide");
p2 = plot(output1.gen_ave, xlims=(0,720), ylims=(0.0,35.0), title="Average Generation", xticks = 0:24:720, rotation=45,legend=:best,label="Average Generation sp1");
plot!(p2,output2.gen_ave, xlims=(0,720), xticks = 0:24:720, rotation=45,legend=:best,label="Average Generation sp2");
p3 = plot(output1.size_ave, xlims=(0,720), ylims=(0.0,3.1), title="Relative Size sp1", xticks = 0:24:720, rotation=45, legend=:topright, label="Relative Size sp1");
p4 = plot(output1.Cq1_ave, xlims=(0,720), ylims=(0.0,1.8e-10), title="Average C&N Quotas (mmol)", xticks = 0:24:720, rotation=45, label = "Cq1 sp1");
plot!(p4,output1.Cq2_ave, xlims=(0,720), xticks = 0:24:720, rotation=45, label="Cq2 sp1");
plot!(p4,output1.Nq_ave, xlims=(0,720), xticks = 0:24:720, rotation=45, label="Nq sp1");
p5 = plot(output2.size_ave, xlims=(0,720), ylims=(0.0,3.1), title="Relative Size sp2", xticks = 0:24:720, rotation=45, legend=:topright, label="Relative Size sp2");
p6= plot(output2.Cq1_ave, xlims=(0,720), ylims=(0.0,1.2e-9), title="Average C&N Quotas (mmol)", xticks = 0:24:720, rotation=45, label = "Cq1 sp2");
plot!(p6,output2.Cq2_ave, xlims=(0,720), xticks = 0:24:720, rotation=45, label="Cq2 sp2");
plot!(p6,output2.Nq_ave, xlims=(0,720), xticks = 0:24:720, rotation=45, label="Nq sp2");
p7 = plot(nutrients.DIN, xlims=(0,720), ylims=(0.0,5.0e-5), title= "DIN&DON (mmol)", xticks = 0:24:720, rotation=45,label="DIN",legend=:topleft);
plot!(p7,nutrients.DON, xlims=(0,720), xticks = 0:24:720, rotation=45, label="DON");
p8 = plot(nutrients.DOC, xlims=(0,720), ylims=(0.0,3.0e-4), title= "POC&DOC (mmol)", xticks = 0:24:720, rotation=45,label="DOC",legend=:topleft);
plot!(p8,nutrients.POC, xlims=(0,720), xticks = 0:24:720, rotation=45, label="POC");
plt=plot(p1,p2,p3,p4,p5,p6,p7,p8,layout=grid(4,2),size=(1200,1000),dpi=200)
savefig(plt,"results/plot.png")




