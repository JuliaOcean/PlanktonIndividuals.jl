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
Compute total gtr (supposed to be 0), and surface vertical tracer flux(supposed to be 0)
Write a brief summary of each time step into a txt file
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
    println(Cio,@sprintf("%3.0f  %.16E  %.16E",t,Σgtrᶜ,TC))
    println(Nio,@sprintf("%3.0f  %.16E  %.16E",t,Σgtrⁿ,TN))
    println(Pio,@sprintf("%3.0f  %.16E  %.16E",t,Σgtrᵖ,TP))
    close(Cio);close(Nio);close(Pio);
end

"""
    write_pop_dynamics(t, agent_num, counts, filepath)
Write a brief summary of population changes at each time step into a txt file
"""
function write_pop_dynamics(t::Int64, pop, counts, filepath)
    POPio = open(filepath*"dynamic_population.txt","a");
    println(POPio,@sprintf("%3.0f  %7.0f  %5.0f  %5.0f  %5.0f",
                           t,pop,counts.divid,counts.graze,counts.death))
    close(POPio);
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
