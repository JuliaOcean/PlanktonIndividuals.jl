"""
    write_nut_nc_each_step(nut, filepath)
Write a NetCDF file of nutrient fields at each time step
Default filepath -> "results/nutrients/nut."*lpad(string(t),4,"0")*".nc"
"""
function write_nut_nc_each_step(nut::NamedTuple, t::Int64, filepath::String)
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    path = filepath*"nutrients/nut."*lpad(t,10,"0")*".nc"
    isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    v1 = defVar(ds, "DIC", nut.DIC.data, ("xC", "yC", "zC"), attrib = C_attr)
    v2 = defVar(ds, "DOC", nut.DOC.data, ("xC", "yC", "zC"), attrib = C_attr)
    v3 = defVar(ds, "POC", nut.POC.data, ("xC", "yC", "zC"), attrib = C_attr)
    v4 = defVar(ds, "NH4", nut.NH4.data, ("xC", "yC", "zC"), attrib = N_attr)
    v5 = defVar(ds, "NO3", nut.NO3.data, ("xC", "yC", "zC"), attrib = N_attr)
    v6 = defVar(ds, "DON", nut.DON.data, ("xC", "yC", "zC"), attrib = N_attr)
    v7 = defVar(ds, "PON", nut.PON.data, ("xC", "yC", "zC"), attrib = N_attr)
    v8 = defVar(ds, "PO4", nut.PO4.data, ("xC", "yC", "zC"), attrib = P_attr)
    v9 = defVar(ds, "DOP", nut.DOP.data, ("xC", "yC", "zC"), attrib = P_attr)
    v10= defVar(ds, "POP", nut.POP.data, ("xC", "yC", "zC"), attrib = P_attr)
    close(ds)
end

"""
    write_nut_nc_alltime(a, DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime)
Write a NetCDF file of nutrient fields for the whole run, especially for 0D configuration
Default filepath -> "results/nutrients.nc"
"""
function write_nut_nc_alltime(DIC, NH4, NO3, PO4, DOC, DON, DOP, POC, PON, POP, nTime,
                              filepath = "./results/nutrients.nc")
    C_attr = Dict("units" => "mmolC/m^3")
    N_attr = Dict("units" => "mmolN/m^3")
    P_attr = Dict("units" => "mmolP/m^3")
    isfile(filepath) && rm(filepath)
    ds = NCDataset(filepath, "c")
    v1 = defVar(ds, "DIC", DIC, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v2 = defVar(ds, "DOC", DOC, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v3 = defVar(ds, "POC", POC, ("xC", "yC", "zC", "T"), attrib = C_attr)
    v4 = defVar(ds, "NH4", NH4, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v5 = defVar(ds, "NO3", NO3, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v6 = defVar(ds, "DON", DON, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v7 = defVar(ds, "PON", PON, ("xC", "yC", "zC", "T"), attrib = N_attr)
    v8 = defVar(ds, "PO4", PO4, ("xC", "yC", "zC", "T"), attrib = P_attr)
    v9 = defVar(ds, "DOP", DOP, ("xC", "yC", "zC", "T"), attrib = P_attr)
    v10= defVar(ds, "POP", POP, ("xC", "yC", "zC", "T"), attrib = P_attr)
    close(ds)
end

"""
    write_nut_cons(g, gtr, nutₜ, vel, t, filepath)
Write a brief summary of nutrients at each time step into a txt file
"""
function write_nut_cons(g::Grids, gtr::NamedTuple, nutₜ::NamedTuple, t::Int64, filepath)
    Σgtrⁿ = sum(interior(gtr.NH4.data, g) .* g.V) +
            sum(interior(gtr.NO3.data, g) .* g.V) +
            sum(interior(gtr.DON.data, g) .* g.V) +
            sum(interior(gtr.PON.data, g) .* g.V)

    Σgtrᶜ = sum(interior(gtr.DIC.data, g) .* g.V) +
            sum(interior(gtr.DOC.data, g) .* g.V) +
            sum(interior(gtr.POC.data, g) .* g.V)

    Σgtrᵖ = sum(interior(gtr.PO4.data, g) .* g.V) +
            sum(interior(gtr.DOP.data, g) .* g.V) +
            sum(interior(gtr.POP.data, g) .* g.V)

    TN = sum(interior(nutₜ.NH4.data, g) .* g.V) +
         sum(interior(nutₜ.NO3.data, g) .* g.V) +
         sum(interior(nutₜ.DON.data, g) .* g.V) +
         sum(interior(nutₜ.PON.data, g) .* g.V)

    TC = sum(interior(nutₜ.DIC.data, g) .* g.V) +
         sum(interior(nutₜ.DOC.data, g) .* g.V) +
         sum(interior(nutₜ.POC.data, g) .* g.V)

    TP = sum(interior(nutₜ.PO4.data, g) .* g.V) +
         sum(interior(nutₜ.DOP.data, g) .* g.V) +
         sum(interior(nutₜ.POP.data, g) .* g.V)

    Cio = open(filepath*"cons_C.txt","a")
    Nio = open(filepath*"cons_N.txt","a")
    Pio = open(filepath*"cons_P.txt","a")

    println(Cio,@sprintf("%4.0f  %.16E  %.16E  %.4f", t, Σgtrᶜ, TC, mean(interior(nutₜ.DOC.data, g))))
    println(Nio,@sprintf("%4.0f  %.16E  %.16E  %.4f  %.4f",
                         t, Σgtrⁿ, TN, mean(interior(nutₜ.NH4.data, g)), mean(interior(nutₜ.NO3.data, g))))
    println(Pio,@sprintf("%4.0f  %.16E  %.16E  %.4f",t, Σgtrᵖ, TP, mean(interior(nutₜ.PO4.data, g))))
    close(Cio);close(Nio);close(Pio);
end
function write_nut_cons(g::Grids, nutₜ::NamedTuple, t::Int64, filepath)
    Cio = open(filepath*"cons_C.txt","a")
    Nio = open(filepath*"cons_N.txt","a")
    Pio = open(filepath*"cons_P.txt","a")
    println(Cio,@sprintf("%4.0f  %.4f",t, mean(interior(nutₜ.DOC.data, g))))
    println(Nio,@sprintf("%4.0f  %.4f  %.4f",t, mean(interior(nutₜ.NH4.data, g)),mean(interior(nutₜ.NO3.data, g))))
    println(Pio,@sprintf("%4.0f  %.4f",t, mean(interior(nutₜ.PO4.data, g))))
    close(Cio);close(Nio);close(Pio);
end

"""
    write_species_dynamics(t, phyts, filepath)
Write a brief summary of each species at each time step into a txt file
"""
function write_species_dynamics(t::Int64, phytos, filepath)
    for i in 1:length(phytos)
        pop = size(phytos[i].data,1)
        gen_ave = mean(phytos[i].data[:,11])
        age_ave = mean(phytos[i].data[:,12])
        size_ave= mean(phytos[i].data[:,5])
        Bm_ave= mean(phytos[i].data[:,6])
        Cq_ave= mean(phytos[i].data[:,7])
        Nq_ave= mean(phytos[i].data[:,8])
        Pq_ave= mean(phytos[i].data[:,9])
        Chl_ave= mean(phytos[i].data[:,10])
        io = open(filepath*"dynamic_species"*lpad(i,3,"0")*".txt","a");
        println(io,@sprintf("%4.0f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E",t,pop,gen_ave,age_ave,size_ave,Bm_ave,Cq_ave,Nq_ave,Pq_ave,Chl_ave))
        close(io);
    end
end

"""
    write_output(individuals,filepath,time)
write model output of individuals at each time step in a binary file
time = model.t
"""
function write_output(phytos::NamedTuple, filepath, time)
    for i in 1:length(phytos)
        path = filepath*"planks/phy"*lpad(time, 10, "0")*"_"*lpad(i,2,"0")*".bin"
        open(path, "w") do io
            serialize(io, phytos[i].data[:,1:12])
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
    mkdir("$res"*"planks/")
    return res
end
