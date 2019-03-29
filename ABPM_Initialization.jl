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

function setup_agents(N::Int64,Cquota::Array,mean::Float64,var::Float64,bdry)
    phyts0 = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
    for i in 1:N
        # agent location
        x = rand(bdry[1,1]*10:bdry[1,2]*10)/10
        y = rand(bdry[2,1]*10:bdry[2,2]*5)/10
        z = rand(bdry[3,1]*10:bdry[3,2]*7.5)/10
        # a normal distribution with mean variance
        radm = max(0.05, rand(Normal(mean,var)))
        gen  = 1
        size = radm
        Cq1  = Cquota[1]*1E-3
        Cq2  = Cquota[1]*radm
        Nq   = 13/106*2*Cq2
        chl  = Cq2*0.4
        sp   = 1
        push!(phyts0,(x=x,y=y,z=z,gen=gen,size=size,Cq1=Cq1,Cq2=Cq2,Nq=Nq,chl=chl,sp=sp))
    end
    for i in N+1:2N
        # agent location
        x = rand(bdry[1,1]*10:bdry[1,2]*10)/10
        y = rand(bdry[2,2]*5:bdry[2,2]*10)/10
        z = rand(bdry[3,1]*10:bdry[3,2]*7.5)/10
        # a normal distribution with mean variance
        radm = max(0.05, rand(Normal(mean,var)))
        gen  = 1
        size = radm
        Cq1  = Cquota[2]*1E-3
        Cq2  = Cquota[2]*radm
        Nq   = 13/106*2*Cq2
        chl  = Cq2*0.4
        sp   = 2
        push!(phyts0,(x=x,y=y,z=z,gen=gen,size=size,Cq1=Cq1,Cq2=Cq2,Nq=Nq,chl=chl,sp=sp))
    end
    B = [phyts0]
    return B
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

function count_num(phyts_a, bdry)
    cells = zeros(bdry[2,2], bdry[1,2], bdry[3,2])
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        x = trunc(Int, phyt.x)
        y = trunc(Int, phyt.y)
        z = trunc(Int, phyt.z)
        cells[y, x, z] = cells[y, x, z] + 1
    end
    return cells
end

function convert_coordinates(phyts, zf, xg, yg, zgrid, xgrid, ygrid)
    for i in 1:size(phyts,1)
	phyt = phyts[i,:]
	z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
	dz = phyt.z - z; dx = phyt.x - x; dy = phyt.y - y;
	phyt.x = xg[x] - dx * xgrid[x];
	phyt.y = yg[y] + dy * ygrid[y];
	phyt.z = zf[z] - dz * zgrid[z];
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
