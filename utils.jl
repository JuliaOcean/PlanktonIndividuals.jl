function read_input(csv::String,myTime::Int64)
    # irradiance(μmol photons/m^2) and temperature
    # will change later to make it able to choose different month
    input = CSV.read(csv)
    days=myTime÷24 
    temp = copy(input.Temp_Aug)
    for i in 1:days-1
        tt = copy(input.Temp_Aug)
        tt = append!(temp,tt)
    end
    IR = copy(input.IR_Aug)
    for i in 1:days-1
        tt = copy(input.IR_Aug)
        tt = append!(IR,tt)
    end
    return temp,IR
end
function read_offline_vels(fieldroot::String)
    u = ncread(fieldroot*"UVEL_big.nc","u"); # zonal current speed, dx, positive to west
    v = ncread(fieldroot*"VVEL_big.nc","v"); # zonal current speed, dx, positive to west
    w = ncread(fieldroot*"WVEL_big.nc","w"); # zonal current speed, dx, positive to west
    vel = velocity_fields(u, v, w)
    return vel
end
function grid_offline(fieldroot::String)
    zF = ncread(fieldroot*"WVEL_big.nc","zF"); # Cell faces depths
    xF = ncread(fieldroot*"UVEL_big.nc","xF"); # Cell faces point long..
    yF = ncread(fieldroot*"VVEL_big.nc","yF"); # Cell faces point lati..
    zC = ncread(fieldroot*"UVEL_big.nc","zC"); # Cell centers depths
    xC = ncread(fieldroot*"WVEL_big.nc","xC"); # Cell centers point long..
    yC = ncread(fieldroot*"WVEL_big.nc","yC"); # Cell centers point lati..
    Nx = length(xF); Ny = length(yF); Nz = length(zF);
    Δz = zF[1:end-1] .- zF[2:end]; # unit: meters
    Δx = (xC[2:end] .- xC[1:end-1]) .* (111.32*cos(π/6)*1000); # unit: meters, at 30N
    Δy = (yC[2:end] .- yC[1:end-1]) .* (111*1000); # unit: meters
    g = grids(xC, yC, zC, xF, yF, zF, Δx, Δy, Δz, Nx, Ny, Nz)
    return g
end

function create_output(B::Array{DataFrame,1})
    output = DataFrame(time=0, gen_ave=mean(B[1].gen), spec_ave = mean(B[1].sp), Cq1_ave=mean(B[1].Cq1), Cq2_ave=mean(B[1].Cq2), Nq_ave=mean(B[1].Nq), size_ave=mean(B[1].size), chl_ave=mean(B[1].chl), Population=size(B[1],1), dvid=0, graz=0)
    return output
end

function write_output(t,CR,output)
    # summary of current step
    gen_ave=mean(CR[1].gen)
    spec_ave=mean(CR[1].sp)
    Cq1_ave=mean(CR[1].Cq1)
    Cq2_ave=mean(CR[1].Cq2)
    Nq_ave=mean(CR[1].Nq)
    size_ave=mean(CR[1].size)
    chl_ave=mean(CR[1].chl)
    push!(output,(time=t, gen_ave=gen_ave, spec_ave=spec_ave, Cq1_ave=Cq1_ave, Cq2_ave=Cq2_ave, Nq_ave=Nq_ave, size_ave=size_ave, chl_ave=chl_ave, Population=size(CR[1],1), dvid=CR[2], graz=CR[3]))
    return output
end

function count_num(phyts_a, grid)
    cells = zeros(grid.Ny, grid.Nx, grid.Nz)
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        x = trunc(Int, phyt.x)
        y = trunc(Int, phyt.y)
        z = trunc(Int, phyt.z)
        cells[y, x, z] = cells[y, x, z] + 1
    end
    return cells
end

function convert_coordinates(phyts, grid)
    for i in 1:size(phyts,1)
    phyt = phyts[i,:]
    z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
    dz = phyt.z - z; dx = phyt.x - x; dy = phyt.y - y;
    phyt.x = grid.xF[x] - dx * grid.Δx[x] ./(1000*111.32*cos(π/6)); # converted to degree
    phyt.y = grid.yF[y] + dy * grid.Δy[y] ./(1000*111.0);
    phyt.z = grid.zF[z] - dz * grid.Δz[z];
    end
end

function sort_species(Bi, B1, B2)
    phyts1 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[],sp=Int64[])
    phyts2 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[],sp=Int64[])
    for j in 1:size(Bi,1)
        if Bi[j,:].sp == 1
            append!(phyts1,Bi[j,:])
        elseif Bi[j,:].sp == 2
            append!(phyts2,Bi[j,:])
        end
    end
    push!(B1,phyts1)
    push!(B2,phyts2)
end
    
function compute_mean_species(B1, B2)
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
    return output1, output2
end
