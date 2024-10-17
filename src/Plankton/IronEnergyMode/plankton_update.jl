function plankton_update!(plank, nuts, rnd, p, plk, diags_spcs, ΔT, t, arch::Architecture, mode::AbstractMode)
    plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch)

    calc_consume!(plk.DIC.data, plk.NH4.data, plk.NO3.data, plk.PO4.data, plk.FeT.data,
                  plank, plank.ac, plank.xi, plank.yi, plank.zi, ΔT, arch)
    ##### diagnostics of processes for each species
    diags_spcs!(diags_spcs, plank, plank.ac, plank.xi, plank.yi, plank.zi, mode, arch)

    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10.0f0,3600.0f0÷ΔT))) == 0.0f0 
        ##### grazing and its diagnostic
        diags_proc!(diags_spcs.graz, plank.graz, plank.ac, plank.xi, plank.yi, plank.zi, arch)

        grazing!(plank, arch, plk, p)

        ###### mortality and its diagnostic
        diags_proc!(diags_spcs.mort, plank.mort, plank.ac, plank.xi, plank.yi, plank.zi, arch)

        mortality!(plank, arch, plk, p)

        ###### cell division diagnostic
        diags_proc!(diags_spcs.dvid, plank.dvid, plank.ac, plank.xi, plank.yi, plank.zi, arch)

        ##### division
        ##### check if the number of individuals exceeded
        ##### do not copy inactive individuals
        plank.dvid .*= plank.ac
        dvidnum = dot(plank.dvid, plank.ac)
        deactive_ind = findall(x -> x == 0.0f0, plank.ac)
        if dvidnum > length(deactive_ind)
            throw(ArgumentError("number of individual exceeds the capacity at timestep $(t/86400.0) days"))
        end
        divide!(plank, nuts, deactive_ind, arch)
        nuts.idc .= 0
    end

    ##### diagnostic for individual distribution
    diags_proc!(diags_spcs.num, plank.ac, plank.ac, plank.xi, plank.yi, plank.zi, arch)
end