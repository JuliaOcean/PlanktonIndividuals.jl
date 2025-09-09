function plankton_update!(phyto, trs, rnd, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)
    plank = phyto.data
    p = phyto.p
    
    ##### grazing and its diagnostic
    diags_proc!(diags_spcs.graz, plank.graz, plank.ac, plank.xi, plank.yi, plank.zi, arch)
    grazing!(plank, arch, plk, p)

    ###### mortality and its diagnostic
    diags_proc!(diags_spcs.mort, plank.mort, plank.ac, plank.xi, plank.yi, plank.zi, arch)
    mortality!(plank, arch, plk, p)

    ##### division
    ##### check if the number of individuals exceeded
    ##### do not copy inactive individuals
    plank.dvid .*= plank.ac
    diags_proc!(diags_spcs.dvid, plank.dvid, plank.ac, plank.xi, plank.yi, plank.zi, arch)
    dvidnum = dot(plank.dvid, plank.ac)
    deactive_ind = findall(x -> x == 0.0f0, plank.ac)
    if dvidnum > length(deactive_ind)
        throw(ArgumentError("number of individual exceeds the capacity at timestep $(t/86400.0) days"))
    end
    divide!(plank, trs, deactive_ind, arch)
    unsafe_free!(deactive_ind)

    ##### phytoplankton physiology processes and interactions with fields
    plankton_growth!(plank, trs, rnd, p, ΔT, t, arch)
    calc_consume!(plk.DIC.data, plank, plank.ac, plank.xi, plank.yi, plank.zi, ΔT, arch)
    
    ##### diagnostics of processes for each species
    diags_spcs!(diags_spcs, phyto, plank.ac, plank.xi, plank.yi, plank.zi, mode, arch)
    diags_proc!(diags_spcs.num, plank.ac, plank.ac, plank.xi, plank.yi, plank.zi, arch)
end