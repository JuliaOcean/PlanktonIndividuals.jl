function plankton_update!(plank, nuts, proc, p, plk, diags_spcs, ΔT, t, arch::Architecture, mode::QuotaMode)
    plankton_growth!(plank, nuts, proc, p, ΔT, t, arch)

    calc_consume!(plk.DIC.data, plk.DOC.data, plk.NH4.data, plk.NO3.data, plk.PO4.data,
                  proc, plank.ac, plank.xi, plank.yi, plank.zi, ΔT, arch)
    ##### diagnostics of processes for each species
    diags_spcs!(diags_spcs, proc, plank, plank.ac, plank.xi, plank.yi, plank.zi, arch)

    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
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
        dvidnum = dot(plank.dvid, plank.ac)
        deactive_ind = findall(x -> x == 0.0, plank.ac)
        if dvidnum > length(deactive_ind)
            throw(ArgumentError("number of individual exceeds the capacity"))
        end
        ##### do not copy inactive individuals
        plank.dvid .*= plank.ac
        divide!(plank, deactive_ind, arch)
    end

    ##### diagnostic for individual distribution
    diags_proc!(diags_spcs.num, plank.ac, plank.ac, plank.xi, plank.yi, plank.zi, arch)
end

function plankton_update!(plank, nuts, proc, p, plk, ::Nothing, ΔT, t, arch::Architecture, mode::QuotaMode)
    plankton_growth!(plank, nuts, proc, p, ΔT, t, arch)

    calc_consume!(plk.DIC.data, plk.DOC.data, plk.NH4.data, plk.NO3.data, plk.PO4.data,
                  proc, plank.ac, plank.xi, plank.yi, plank.zi, ΔT, arch)

    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
        ##### grazing
        grazing!(plank, arch, plk, p)

        ###### mortality
        mortality!(plank, arch, plk, p)

        ##### division
        ##### check if the number of individuals exceeded
        dvidnum = dot(plank.dvid, plank.ac)
        deactive_ind = findall(x -> x == 0.0, plank.ac)
        if dvidnum > length(deactive_ind)
            throw(ArgumentError("number of individual exceeds the capacity"))
        end
        ##### do not copy inactive individuals
        plank.dvid .*= plank.ac
        divide!(plank, deactive_ind, arch)
    end
end