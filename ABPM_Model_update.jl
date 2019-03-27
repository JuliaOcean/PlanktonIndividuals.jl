using DataFrames, CSV, LaTeXStrings, Plots, NetCDF
using Random
using Distributions
cd("/home/zhenwu/ABPM_3D")
include("ABPM_Parameters.jl")
include("ABPM_Functions.jl")
include("ABPM_Initialization.jl")
# Read input files
myTime = 720 # runtime, time step: 1h
temp,IR = read_input("T_IR.csv",myTime);
fieldroot = "/nobackup1b/users/zhenwu/vel_fields/";
u = ncread(fieldroot*"UVEL_Mar.nc","u"); # zonal current speed, dx, positive to west
v = ncread(fieldroot*"VVEL_Mar.nc","v"); # zonal current speed, dx, positive to west
w = ncread(fieldroot*"WVEL_Mar.nc","w"); # zonal current speed, dx, positive to west
zf = ncread(fieldroot*"WVEL_Mar.nc","zC"); # Cell faces depths
xg = ncread(fieldroot*"WVEL_Mar.nc","xC"); # Cell corner point long..
yg = ncread(fieldroot*"WVEL_Mar.nc","yC"); # Cell corner point lati..
zgrid = zf[1:end-1] .- zf[2:end];
xgrid = xg[2:end] .- xg[1:end-1];
ygrid = yg[2:end] .- yg[1:end-1];
# SETUP
bdry = [[1 400]; [1 300]; [1 30]]; # boundaries of grids
N = 40000   # Number of initial individuals of each species
Nsp = 2     # Number of species
B=setup_agents(N,Cquota,1.1,0.18,bdry) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nutrients = DataFrame(DIN=5.0e-6, DOC=0.0, DON=6.0e-5, POC=0.0, PON=0.0);# μmol
# model update
for t in 1:myTime
    phyts_a = copy(B[t]) # read data from last time step
    agent_move(phyts_a,bdry,u,v,w,xgrid,ygrid,zgrid,t) # agents advection and  convection
    cell_num = count_num(phyts_a, bdry)
    CR=update(t, phyts_a, nutrients, IR, temp, cell_num) # model update, return value: phyts_b, dvid_ct, and graz_ct
    push!(B,CR[1])
    write_output(t,CR,output)
    println(output[t,:]) # save current output
    println(nutrients[t,:])
end
B1 = []; B2 = [];
for i in 1:720
    phyts1 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[],sp=Int64[])
    phyts2 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[],sp=Int64[])
    for j in 1:size(B[i],1)
        if B[i][j,:].sp == 1
            append!(phyts1,B[i][j,:])
        elseif B[i][j,:].sp == 2
            append!(phyts2,B[i][j,:])
        end
    end
    push!(B1,phyts1)
    push!(B2,phyts2)
end
output1 = DataFrame(time=Int64[], gen_ave=Float64[], Cq1_ave=Float64[], Cq2_ave=Float64[], Nq_ave=Float64[], size_ave=Float64[], chl_ave=Float64[], Population=Int64[]);
output2 = DataFrame(time=Int64[], gen_ave=Float64[], Cq1_ave=Float64[], Cq2_ave=Float64[], Nq_ave=Float64[], size_ave=Float64[], chl_ave=Float64[], Population=Int64[]);
for i in 1:720
    gen_ave1=mean(B1[i].gen)
    Cq1_ave1=mean(B1[i].Cq1)
    Cq2_ave1=mean(B1[i].Cq2)
    Nq_ave1=mean(B1[i].Nq)
    size_ave1=mean(B1[i].size)
    chl_ave1=mean(B1[i].chl)
    push!(output1,(time=i, gen_ave=gen_ave1, Cq1_ave=Cq1_ave1, Cq2_ave=Cq2_ave1, Nq_ave=Nq_ave1, size_ave=size_ave1, chl_ave=chl_ave1, Population=size(B1[i],1)))
    gen_ave2=mean(B2[i].gen)
    Cq1_ave2=mean(B2[i].Cq1)
    Cq2_ave2=mean(B2[i].Cq2)
    Nq_ave2=mean(B2[i].Nq)
    size_ave2=mean(B2[i].size)
    chl_ave2=mean(B2[i].chl)
    push!(output2,(time=i, gen_ave=gen_ave2, Cq1_ave=Cq1_ave2, Cq2_ave=Cq2_ave2, Nq_ave=Nq_ave2, size_ave=size_ave2, chl_ave=chl_ave2, Population=size(B2[i],1)))
end
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
 close(f)

