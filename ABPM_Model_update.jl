using DataFrames, CSV, LaTeXStrings, Plots, NetCDF, ProgressMeter
using Random
using Distributions
cd("/home/zhenwu/ABPM_3D")
include("ABPM_Parameters.jl")
include("ABPM_Functions.jl")
include("ABPM_Initialization.jl")
# Read input files
nTime = 720 # number of time steps
deltaT = 3600 # time step : 1h
temp,IR = read_input("T_IR.csv",nTime);
fieldroot = "/nobackup1b/users/zhenwu/vel_fields/";
#u = ncread(fieldroot*"UVEL_Mar.nc","u"); # zonal current speed, dx, positive to west
#v = ncread(fieldroot*"VVEL_Mar.nc","v"); # zonal current speed, dx, positive to west
#w = ncread(fieldroot*"WVEL_Mar.nc","w"); # zonal current speed, dx, positive to west
#zf = ncread(fieldroot*"WVEL_Mar.nc","zC"); # Cell faces depths
#xg = ncread(fieldroot*"WVEL_Mar.nc","xC"); # Cell corner point long..
#yg = ncread(fieldroot*"WVEL_Mar.nc","yC"); # Cell corner point lati..
vel = read_offline_vels(fieldroot);
g = grid_offline(fieldroot);
# SETUP
N = 60000   # Number of initial individuals of each species
Nsp = 2     # Number of species
B=setup_agents(N,Cquota,1.1,0.18,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nutrients = DataFrame(DIN=5.0e-6, DOC=0.0, DON=6.0e-5, POC=0.0, PON=0.0);# μmol
# model update
for t in 1:nTime
    phyts_a = copy(B[t]) # read data from last time step
    agent_move(phyts_a,vel,g,t,deltaT) 
    cell_num = count_num(phyts_a, g)
    CR=update(t, deltaT, phyts_a, nutrients, IR, temp, cell_num) # model update, return value: phyts_b, dvid_ct, and graz_ct
    push!(B,CR[1])
    write_output(t,CR,output)
    println(output[t,:]) # save current output
    println(nutrients[t,:])
    convert_coordinates(B[t],g) # convert grids to lon, lat and depth
end
B1 = []; B2 = [];
@showprogress 1 "Computing..." for i in 1:720
    sort_species(B[i], B1, B2)
end
output1, output2 = compute_mean_species(B1, B2)
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
 JLD.@write f nutrients
 JLD.@write f bdry
 JLD.@write f zf
 JLD.@write f xg
 JLD.@write f yg
 close(f)

