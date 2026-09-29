function colony_update!(colony, trs, rnd, plk, diags_colony, ΔT, t, arch::Architecture)
    ##### population dynamics
    ##### species 1
    sp1 = colony.spcs.sp1.data
    sp1_p = colony.spcs.sp1.p
    diags_proc!(diags_colony.sp1.graz, sp1.graz, sp1.ac, sp1.xi, sp1.yi, sp1.zi, arch)
    grazing!(sp1, arch, plk, sp1_p)
        
    diags_proc!(diags_colony.sp1.mort, sp1.mort, sp1.ac, sp1.xi, sp1.yi, sp1.zi, arch)
    mortality!(sp1, arch, plk, sp1_p)
        
    sp1.dvid .*= sp1.ac
    diags_proc!(diags_colony.sp1.dvid, sp1.dvid, sp1.ac, sp1.xi, sp1.yi, sp1.zi,arch)
        
    dvidnum = dot(sp1.dvid, sp1.ac)
    deactive_ind = findall(isequal(false), sp1.ac)
    if dvidnum > length(deactive_ind)
        throw(ArgumentError("number of individual exceeds the capacity at timestep $(t/86400.0) days"))
    end
    divide!(sp1, trs, deactive_ind, arch)
    unsafe_free!(deactive_ind)

    diags_proc!(diags_colony.sp1.num, sp1.ac, sp1.ac, sp1.xi, sp1.yi, sp1.zi, arch)

    for sp in eachindex(colony.spcs)[2:end]
        plank = colony.spcs[sp].data
        p = colony.spcs[sp].p

        grazing!(plank, arch, plk, p)
        mortality!(plank, arch, plk, p)
        
        plank.dvid .*= plank.ac
        dvidnum = dot(plank.dvid, plank.ac)
        deactive_ind = findall(isequal(false), plank.ac)
        if dvidnum > length(deactive_ind)
            throw(ArgumentError("number of individual exceeds the capacity at timestep $(t/86400.0) days"))
        end
        divide!(plank, trs, deactive_ind, arch)
        unsafe_free!(deactive_ind)
    end
    
    ##### interaction between different cells
    for pair in colony.intac
        sp1 = colony.spcs[pair[1]].data
        sp1_p = colony.spcs[pair[1]].p
        sp2 = colony.spcs[pair[2]].data
        calc_material_exchange!(sp1, sp2, sp1_p, arch)
        update_material_exchange!(sp1, sp2, ΔT, arch)
    end
    
    ##### phytoplankton physiological processes and interactions with fields
    for sp in eachindex(colony.spcs)
        plank = colony.spcs[sp].data
        p = colony.spcs[sp].p
        colony_plankton_growth!(plank, trs, p, ΔT, arch)
        calc_consume!(plk.DIC.data, plk.NH4.data, plk.NO3.data, plk.PO4.data, plk.DFe.data, plk.O2.data,
                      plank, plank.ac, plank.xi, plank.yi, plank.zi, ΔT, arch)
    end # physiological processes

    calc_plankton_population_dynamics(sp1, trs, rnd, sp1_p, ΔT, t, arch)
    for sp in eachindex(colony.spcs)[2:end]
        plank = colony.spcs[sp].data
        plank.dvid .= sp1.dvid
        plank.graz .= sp1.graz
        plank.mort .= sp1.mort
    end    
  
    ##### diagnostics of physiological processes for each species
    diags_colony!(diags_colony, colony.spcs, arch)
end