"""
    write_nut_nc_each_step(g, nut, t)
Write a NetCDF file of nutrient fields at each time step
Default filepath -> "results/nutrients/nut."*lpad(string(t),4,"0")*".nc"
"""
function write_nut_nc_each_step(g::grids, nut::nutrient_fields, t::Int64, filepath::String)
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "NH4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "NO3", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "PO4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "DOP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=N_attr);
    nccreate(filepath, "POP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, atts=P_attr);
    ncwrite(nut.DIC,filepath,"DIC"); ncwrite(nut.NH4,filepath,"NH4");
    ncwrite(nut.NO3,filepath,"NO3"); ncwrite(nut.PO4,filepath,"PO4");
    ncwrite(nut.DOC,filepath,"DOC"); ncwrite(nut.DON,filepath,"DON"); ncwrite(nut.DOP,filepath,"DOP");
    ncwrite(nut.POC,filepath,"POC"); ncwrite(nut.PON,filepath,"PON"); ncwrite(nut.POP,filepath,"POP");
    nothing
end

"""
    write_nut_nc_alltime(a, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime)
Write a NetCDF file of nutrient fields for the whole run, especially for 0D configuration
Default filepath -> "results/nutrients.nc"
"""
function write_nut_nc_alltime(g::grids, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime,
                              filepath = "./results/nutrients.nc")
    tt = collect(1:nTime);
    xC_attr = Dict("longname" => "Locations of the cell centers in the x-direction.", "units" => "m")
    yC_attr = Dict("longname" => "Locations of the cell centers in the y-direction.", "units" => "m")
    zC_attr = Dict("longname" => "Locations of the cell centers in the z-direction.", "units" => "m")
    T_attr = Dict("longname" => "Time", "units" => "H")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    nccreate(filepath, "DIC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "NH4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "NO3", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "PO4", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    nccreate(filepath, "DOC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "DON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "DOP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    nccreate(filepath, "POC", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=C_attr);
    nccreate(filepath, "PON", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=N_attr);
    nccreate(filepath, "POP", "xC", g.xC[:,1], xC_attr, "yC", g.yC[1,:], yC_attr, "zC", g.zC, zC_attr, "T", tt, T_attr, atts=P_attr);
    ncwrite(DIC,filepath,"DIC"); ncwrite(NH4,filepath,"NH4");
    ncwrite(NO3,filepath,"NO3"); ncwrite(PO4,filepath,"PO4");
    ncwrite(DOC,filepath,"DOC"); ncwrite(DON,filepath,"DON"); ncwrite(DOP,filepath,"DOP");
    ncwrite(POC,filepath,"POC"); ncwrite(PON,filepath,"PON"); ncwrite(POP,filepath,"POP");
    ncclose(filepath)
    return nothing
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
function write_pop_dynamics(t::Int64, phyts, counts, filepath)
    pop = size(phyts,2)
    gen_ave = mean(phyts[11,:])
    POPio = open(filepath*"dynamic_population.txt","a");
    println(POPio,@sprintf("%4.0f  %6.0f  %1.2f  %4.0f  %4.0f  %4.0f",t,pop,gen_ave,counts.divid,counts.graze,counts.death))
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
        phyt_sp[i] = reshape(phyt_sp[i],size(phyts,1),Int(length(phyt_sp[i])/size(phyts,1)))
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
function write_output(individuals, filepath, time)
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
