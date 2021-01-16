"""
    write_nut_nc_each_step(nut, t, filepath)
Write a NetCDF file of nutrient fields at each time step

Keyword Arguments
=================
- `nut`: `NamedTuple` of nutrient tracers.
- `t`: Current time of `model` in second, usually starting from 0.
- `filepath`: The file path to store NetCDF files.
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

##### write a brief summary of nutrients at each time step into a txt file
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

    day = t÷86400
    hour = t%86400/3600

    Cio = open(filepath*"cons_C.txt","a")
    Nio = open(filepath*"cons_N.txt","a")
    Pio = open(filepath*"cons_P.txt","a")

    println(Cio,@sprintf("%3.0f  %2.2f  %.16E  %.16E  %.4f", day, hour, Σgtrᶜ, TC, mean(interior(nutₜ.DOC.data, g))))
    println(Nio,@sprintf("%3.0f  %2.2f  %.16E  %.16E  %.4f  %.4f",
                         day, hour, Σgtrⁿ, TN, mean(interior(nutₜ.NH4.data, g)), mean(interior(nutₜ.NO3.data, g))))
    println(Pio,@sprintf("%3.0f  %2.2f  %.16E  %.16E  %.4f",day, hour, Σgtrᵖ, TP, mean(interior(nutₜ.PO4.data, g))))
    close(Cio);close(Nio);close(Pio);
end
function write_nut_cons(g::Grids, nutₜ::NamedTuple, t::Int64, filepath)
    Cio = open(filepath*"cons_C.txt","a")
    Nio = open(filepath*"cons_N.txt","a")
    Pio = open(filepath*"cons_P.txt","a")
    day = t÷86400
    hour = t%86400/3600
    println(Cio,@sprintf("%3.0f  %2.2f  %.4f",day, hour, mean(interior(nutₜ.DOC.data, g))))
    println(Nio,@sprintf("%3.0f  %2.2f  %.4f  %.4f",day, hour, mean(interior(nutₜ.NH4.data, g)),mean(interior(nutₜ.NO3.data, g))))
    println(Pio,@sprintf("%3.0f  %2.2f  %.4f",day, hour, mean(interior(nutₜ.PO4.data, g))))
    close(Cio);close(Nio);close(Pio);
end

##### write a brief summary of each species at each time step into a txt file
function write_species_dynamics(t::Int64, phytos, filepath)
    for i in 1:length(phytos)
        pop = dot(phytos[i].data.ac, phytos[i].data.ac)
        gen_ave =  dot(phytos[i].data.gen, phytos[i].data.ac) / pop
        age_ave =  dot(phytos[i].data.age, phytos[i].data.ac) / pop
        size_ave=  dot(phytos[i].data.Sz,  phytos[i].data.ac) / pop
        Bm_ave  =  dot(phytos[i].data.Bm,  phytos[i].data.ac) / pop
        Cq_ave  =  dot(phytos[i].data.Cq,  phytos[i].data.ac) / pop
        Nq_ave  =  dot(phytos[i].data.Nq,  phytos[i].data.ac) / pop
        Pq_ave  =  dot(phytos[i].data.Pq,  phytos[i].data.ac) / pop
        Chl_ave =  dot(phytos[i].data.chl, phytos[i].data.ac) / pop
        day = t÷86400
        hour = t%86400/3600
        io = open(filepath*"dynamic_species"*lpad(i,3,"0")*".txt","a");
        println(io,@sprintf("%3.0f  %2.2f  %6.0f  %1.2f  %1.2f  %1.2f  %.8E  %.8E  %.8E  %.8E  %.8E",
                            day,hour,pop,gen_ave,age_ave,size_ave,Bm_ave,Cq_ave,Nq_ave,Pq_ave,Chl_ave))
        close(io);
    end
end

"""
    write_individuals_to_bin(phytos, filepath, t)
write model output of individuals at each time step to a binary file

Keyword Arguments
=================
- `phytos`: `NamedTuple` of a list of `individual` species.
- `filepath`: The file path to store JLD2 files.
- `t`: Current time of `model` in second, usually starting from 0.
- `atts` (optional): attributes of individuals to save, default `(:x, :y, :z)`
"""
function write_individuals_to_jld2(phytos::NamedTuple, filepath, t; 
                                   atts = (:x, :y, :z))
    jldopen(filepath*"individuals.jld2", "a+") do file
        for sp in keys(phytos)
            spi = NamedTuple{atts}([getproperty(phytos[sp].data, att) for att in atts])
            for att in atts
                file[lpad(t, 8, "0")*"/"*string(sp)*"/"*string(att)] = Array(spi[att])
            end
        end
    end
end

"""
    write_diags_to_jld2(diags, filepath, t, ncounts)
write model output of individuals at each time step to a binary file

Keyword Arguments
=================
- `diags`: `NamedTuple` of a list of diagnostics at current time step.
- `filepath`: The file path to store JLD2 files.
- `t`: Current time of `model` in second, usually starting from 0.
- `ncounts`: the number of time steps included in each diagnostic
"""
function write_diags_to_jld2(diags, filepath, time, ncounts)
    t = time ÷ 3600 + 1 # time in hour from start
    jldopen(filepath*"diags.jld2", "a+") do file
        for key in keys(diags.tr)
            file[lpad(t, 8, "0")*"/nut/"*string(key)] = Array(diags.tr[key]) ./ ncounts
        end
        for sp in keys(diags.spcs)
            for proc in keys(diags.spcs[sp])
                file[lpad(t, 8, "0")*"/"*string(sp)*"/"*string(proc)] = Array(diags.spcs[sp][proc]) ./ ncounts 
            end
        end
    end
    ##### zeros diags
    for tr in diags.tr
        tr .= 0.0
    end
    for sp in diags.spcs
        for proc in sp
            proc .= 0.0
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
