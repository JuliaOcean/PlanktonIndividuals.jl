"""
    write_nut_nc_each_step(g, nut, t)
Write a NetCDF file of nutrient fields at each time step
Default filepath -> "results/nutrients/nut."*lpad(string(t),4,"0")*".nc"
"""
function write_nut_nc_each_step(g::grids, nut::nutrient_fields, t::Int64, filepath::String)
    xC = g.xC[2:end-1]; yC = g.yC[2:end-1]; zC = g.zC[2:end-1]
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    ds = NCDataset(filepath, "c")
    defDim(ds, "xC", size(xC,1))
    defDim(ds, "yC", size(yC,1))
    defDim(ds, "zC", size(zC,1))
    v1 = defVar(ds, "DIC", Float64, ("xC", "yC", "zC"), attrib = C_attr)
    v2 = defVar(ds, "DOC", Float64, ("xC", "yC", "zC"), attrib = C_attr)
    v3 = defVar(ds, "POC", Float64, ("xC", "yC", "zC"), attrib = C_attr)
    v4 = defVar(ds, "NH4", Float64, ("xC", "yC", "zC"), attrib = N_attr)
    v5 = defVar(ds, "NO3", Float64, ("xC", "yC", "zC"), attrib = N_attr)
    v6 = defVar(ds, "DON", Float64, ("xC", "yC", "zC"), attrib = N_attr)
    v7 = defVar(ds, "PON", Float64, ("xC", "yC", "zC"), attrib = N_attr)
    v8 = defVar(ds, "PO4", Float64, ("xC", "yC", "zC"), attrib = P_attr)
    v9 = defVar(ds, "DOP", Float64, ("xC", "yC", "zC"), attrib = P_attr)
    v10= defVar(ds, "POP", Float64, ("xC", "yC", "zC"), attrib = P_attr)

    v1[:,:,:] = nut.DIC; v2[:,:,:] = nut.DOC; v3[:,:,:] = nut.POC;
    v4[:,:,:] = nut.NH4; v5[:,:,:] = nut.NO3; v6[:,:,:] = nut.DON; v7[:,:,:] = nut.PON;
    v8[:,:,:] = nut.PO4; v9[:,:,:] = nut.DOP; v10[:,:,:] = nut.POP;

    close(ds)
end

"""
    write_nut_nc_alltime(a, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime)
Write a NetCDF file of nutrient fields for the whole run, especially for 0D configuration
Default filepath -> "results/nutrients.nc"
"""
function write_nut_nc_alltime(g::grids, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime,
                              filepath = "./results/nutrients.nc")
    xC = g.xC[2:end-1]; yC = g.yC[2:end-1]; zC = g.zC[2:end-1]
    tt = collect(1:nTime);
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    ds = NCDataset(filepath, "c")
    defDim(ds, "xC", size(xC,1))
    defDim(ds, "yC", size(yC,1))
    defDim(ds, "zC", size(zC,1))
    defDim(ds, "T", size(tt,1))
    v1 = defVar(ds, "DIC", Float64, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v2 = defVar(ds, "DOC", Float64, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v3 = defVar(ds, "POC", Float64, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v4 = defVar(ds, "NH4", Float64, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v5 = defVar(ds, "NO3", Float64, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v6 = defVar(ds, "DON", Float64, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v7 = defVar(ds, "PON", Float64, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v8 = defVar(ds, "PO4", Float64, ("xC", "yC", "zC", "T"), attrib = P_attr)
    v9 = defVar(ds, "DOP", Float64, ("xC", "yC", "zC", "T"), attrib = P_attr)
    v10= defVar(ds, "POP", Float64, ("xC", "yC", "zC", "T"), attrib = P_attr)

    v1[:,:,:,:] = nut.DIC; v2[:,:,:,:] = nut.DOC; v3[:,:,:,:] = nut.POC;
    v4[:,:,:,:] = nut.NH4; v5[:,:,:,:] = nut.NO3; v6[:,:,:,:] = nut.DON; v7[:,:,:,:] = nut.PON;
    v8[:,:,:,:] = nut.PO4; v9[:,:,:,:] = nut.DOP; v10[:,:,:,:] = nut.POP;

    close(ds)
end

"""
    write_nut_cons(g, gtr, nutₜ, vel, t, filepath)
Write a brief summary of nutrients at each time step into a txt file
"""
function write_nut_cons(g::grids, gtr::nutrient_fields, nutₜ::nutrient_fields, t::Int64, filepath)
    Σgtrⁿ = sum(gtr.NH4 .* g.V)+sum(gtr.NO3 .* g.V)+sum(gtr.DON .* g.V)+sum(gtr.PON .* g.V)
    Σgtrᶜ = sum(gtr.DIC .* g.V)+sum(gtr.DOC .* g.V)+sum(gtr.POC .* g.V)
    Σgtrᵖ = sum(gtr.PO4 .* g.V)+sum(gtr.DOP .* g.V)+sum(gtr.POP .* g.V)
    TC = sum((nutₜ.DIC .+ nutₜ.DOC .+ nutₜ.POC) .* g.V)
    TN = sum((nutₜ.NH4 .+ nutₜ.NO3 .+ nutₜ.DON .+ nutₜ.PON) .* g.V)
    TP = sum((nutₜ.PO4 .+ nutₜ.DOP .+ nutₜ.POP) .* g.V)
    Cio = open(filepath*"cons_C.txt","a"); Nio = open(filepath*"cons_N.txt","a");
    Pio = open(filepath*"cons_P.txt","a");
    println(Cio,@sprintf("%4.0f  %.16E  %.16E  %.4f",t,Σgtrᶜ,TC,mean(nutₜ.DOC)))
    println(Nio,@sprintf("%4.0f  %.16E  %.16E  %.4f  %.4f",
                         t,Σgtrⁿ,TN,mean(nutₜ.NH4),mean(nutₜ.NO3)))
    println(Pio,@sprintf("%4.0f  %.16E  %.16E  %.4f",t,Σgtrᵖ,TP,mean(nutₜ.PO4)))
    close(Cio);close(Nio);close(Pio);
end
function write_nut_cons(g::grids, nutₜ::nutrient_fields, t::Int64, filepath)
    TC = sum((nutₜ.DIC .+ nutₜ.DOC .+ nutₜ.POC) .* g.V)
    TN = sum((nutₜ.NH4 .+ nutₜ.NO3 .+ nutₜ.DON .+ nutₜ.PON) .* g.V)
    TP = sum((nutₜ.PO4 .+ nutₜ.DOP .+ nutₜ.POP) .* g.V)
    Cio = open(filepath*"cons_C.txt","a"); Nio = open(filepath*"cons_N.txt","a");
    Pio = open(filepath*"cons_P.txt","a");
    println(Cio,@sprintf("%4.0f  %.16E  %.4f",t,TC,mean(nutₜ.DOC)))
    println(Nio,@sprintf("%4.0f  %.16E  %.4f  %.4f",t,TN,mean(nutₜ.NH4),mean(nutₜ.NO3)))
    println(Pio,@sprintf("%4.0f  %.16E  %.4f",t,TP,mean(nutₜ.PO4)))
    close(Cio);close(Nio);close(Pio);
end

"""
    write_pop_dynamics(t, phyts, counts, filepath)
Write a brief summary of population changes at each time step into a txt file
"""
function write_pop_dynamics(t::Int64, counts, filepath)
    POPio = open(filepath*"dynamic_population.txt","a");
    for i in 1:length(counts.divid)
        println(POPio,@sprintf("%4.0f  %4.0f  %4.0f  %4.0f",t,counts.divid[i],counts.graze[i],counts.death[i]))
    end
    close(POPio);
end

"""
    sort_species(phyts, Nsp)
separate different species in different arrays
"""
function sort_species(phyts,Nsp)
    phyt_sp=[]
    for i in 1:Nsp
        push!(phyt_sp,Real[])
    end
    for i in 1:size(phyts,2)
        phyt = phyts[:,i]
        sp = Int(phyt[10])
        append!(phyt_sp[sp],phyt)
    end
    for i in 1:Nsp
        phyt_sp[i] = reshape(phyt_sp[i],size(phyts,1),Int(length(phyt_sp[i])/size(phyts,1)))
    end
    return phyt_sp
end

"""
    write_species_dynamics(t, phyts, filepath)
Write a brief summary of each species at each time step into a txt file
"""
function write_species_dynamics(t::Int64, phyt_sp, filepath)
    for i in 1:size(phyt_sp,1)
        pop = size(phyt_sp[i],2)
        gen_ave = mean(phyt_sp[i][11,:])
        age_ave = mean(phyt_sp[i][12,:])
        size_ave= mean(phyt_sp[i][4,:])
        Bm_ave= mean(phyt_sp[i][5,:])
        Cq_ave= mean(phyt_sp[i][6,:])
        Nq_ave= mean(phyt_sp[i][7,:])
        Pq_ave= mean(phyt_sp[i][8,:])
        Chl_ave= mean(phyt_sp[i][9,:])
        io = open(filepath*"dynamic_species"*lpad(i,3,"0")*".txt","a");
        println(io,@sprintf("%4.0f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E",t,pop,gen_ave,age_ave,size_ave,Bm_ave,Cq_ave,Nq_ave,Pq_ave,Chl_ave))
        close(io);
    end
end

"""
    write_output(individuals,filepath,time)
write model output of individuals at each time step in a binary file
time = model.t*ΔT
"""
function write_output(individuals::individuals, filepath, time)
    phytos = individuals.phytos
    path = filepath*"phy_"*lpad(time, 10, "0")*".bin"
    if individuals.zoos == nothing
        open(path, "w") do io
            serialize(io, phytos)
        end
    else
        open(path, "w") do io
            serialize(io, phytos)
        end
        path_zoo = filepath*"zoo_"*lpad(time, 10, "0")*".bin"
        open(path_zoo, "w") do io
            serialize(io, individuals.zoos)
        end
    end
end
function write_output(phyts_sp::Array, filepath, time)
    for i in 1:size(phyts_sp,1)
        path = filepath*"phy"*lpad(time, 10, "0")*"_"*lpad(i,2,"0")*".bin"
        open(path, "w") do io
            serialize(io, phyts_sp[i])
        end
    end
end

"""
    PrepRunDir(res::String="results/")
Create `res/` folder if needed. Remove old files from it if needed.
"""
function PrepRunDir(res::String="./results/")
    isdir(res) && rm(res, recursive=true)
    mkdir(res)
    mkdir("$res"*"nutrients/")
    return res
end
