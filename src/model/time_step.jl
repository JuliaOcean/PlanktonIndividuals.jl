"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_Struct, ΔT, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing

    ##### plankton advection
    for sp in keys(model.individuals.phytos)
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd, model.params["κhP"], ΔT)
        periodic_domain!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid)

        plankton_advectionRK4!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)

        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)

        ##### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.grid, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        gen_rand_plk!(model.timestepper.rnd, model.arch)
        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts, 
                         model.individuals.phytos[sp].proc, model.timestepper.rnd, 
                         model.timestepper.par, model.timestepper.pop, model.input.temp[:,:,:,clock], 
                         model.nutrients, model.grid, model.individuals.phytos[sp].p, ΔT, model.t)

        calc_consume!(model.timestepper.plk.DIC.data, model.timestepper.plk.DOC.data, 
                      model.timestepper.plk.NH4.data, model.timestepper.plk.NO3.data, 
                      model.timestepper.plk.PO4.data, model.individuals.phytos[sp].proc, 
                      model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      ΔT, model.grid, model.arch)
        ##### diagnostics of processes for each species
        diags_spcs!(model.diags.spcs[sp], model.individuals.phytos[sp].proc, model.individuals.phytos[sp].data,
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                    model.grid, model.arch)

        ##### grazing and its diagnostic
        diags_graz!(model.diags.spcs[sp].graz, model.individuals.phytos[sp].data.graz,
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                    model.grid, model.arch)
        zero_tmp!(model.timestepper.tmp)
        grazing!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.arch, model.grid, model.timestepper.plk, model.individuals.phytos[sp].p)

        ###### mortality and its diagnostic
        diags_mort!(model.diags.spcs[sp].mort, model.individuals.phytos[sp].data.mort,
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                    model.grid, model.arch)

        zero_tmp!(model.timestepper.tmp)
        mortality!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.arch, model.grid, model.timestepper.plk, model.individuals.phytos[sp].p)

        ###### cell division diagnostic
        diags_dvid!(model.diags.spcs[sp].dvid, model.individuals.phytos[sp].data.dvid,
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                    model.grid, model.arch)

        ##### division
        zero_tmp!(model.timestepper.tmp)
        dvidnum = dot(model.individuals.phytos[sp].data.dvid, model.individuals.phytos[sp].data.ac)
        divide!(model.individuals.phytos[sp].data, model.timestepper.tmp, dvidnum, model.arch)

        ##### tidy up model.individuals.phytos[sp].data, copy to tmp, to the end of divided individuals
        get_tind!(model.individuals.phytos[sp].data.idx, model.individuals.phytos[sp].data.ac)
        model.individuals.phytos[sp].data.idx .= model.individuals.phytos[sp].data.idx .+ dvidnum*2
        copyto_tmp!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.individuals.phytos[sp].data.ac, Int.(model.individuals.phytos[sp].data.idx), 
                    false, model.arch)

        ##### copy back to model.individuals.phytos[sp].data
        zero_tmp!(model.individuals.phytos[sp].data)
        get_tind!(model.timestepper.tmp.idx, model.timestepper.tmp.ac)
        copyto_tmp!(model.timestepper.tmp, model.individuals.phytos[sp].data, 
                    model.timestepper.tmp.ac, Int.(model.timestepper.tmp.idx), 
                    false, model.arch)

        ##### diagnostic for individual distribution
        diags_num!(model.diags.spcs[sp].num, model.individuals.phytos[sp].data.ac, 
                   Int.(model.individuals.phytos[sp].data.xi), Int.(model.individuals.phytos[sp].data.yi), 
                   Int.(model.individuals.phytos[sp].data.zi), model.grid, model.arch)

    end
    write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    @inbounds model.diags.tr.PAR .+= model.timestepper.par
    for key in keys(model.diags.tr)
        if key in keys(model.nutrients)
            @inbounds model.diags.tr[key] .+= model.nutrients[key].data
        end
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end

function PI_TimeStep!(model::Model_Struct, ΔT)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing

    ##### plankton advection
    for sp in keys(model.individuals.phytos)
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd, model.params["κhP"], ΔT)
        periodic_domain!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid)

        plankton_advectionRK4!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)

        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)

        ##### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.grid, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        gen_rand_plk!(model.timestepper.rnd, model.arch)
        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts, 
                         model.individuals.phytos[sp].proc, model.timestepper.rnd, 
                         model.timestepper.par, model.timestepper.pop, model.input.temp[:,:,:,clock], 
                         model.nutrients, model.grid, model.individuals.phytos[sp].p, ΔT, model.t)

        calc_consume!(model.timestepper.plk.DIC.data, model.timestepper.plk.DOC.data, 
                      model.timestepper.plk.NH4.data, model.timestepper.plk.NO3.data, 
                      model.timestepper.plk.PO4.data, model.individuals.phytos[sp].proc, 
                      model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      ΔT, model.grid, model.arch)

        ##### grazing
        zero_tmp!(model.timestepper.tmp)
        grazing!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.arch, model.grid, model.timestepper.plk, model.individuals.phytos[sp].p)

        ###### mortality
        zero_tmp!(model.timestepper.tmp)
        mortality!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.arch, model.grid, model.timestepper.plk, model.individuals.phytos[sp].p)

        ###### cell division
        zero_tmp!(model.timestepper.tmp)
        dvidnum = dot(model.individuals.phytos[sp].data.dvid, model.individuals.phytos[sp].data.ac)
        divide!(model.individuals.phytos[sp].data, model.timestepper.tmp, dvidnum, model.arch)

        ##### tidy up model.individuals.phytos[sp].data, copy to tmp, to the end of divided individuals
        get_tind!(model.individuals.phytos[sp].data.idx, model.individuals.phytos[sp].data.ac)
        model.individuals.phytos[sp].data.idx .= model.individuals.phytos[sp].data.idx .+ dvidnum*2
        copyto_tmp!(model.individuals.phytos[sp].data, model.timestepper.tmp, 
                    model.individuals.phytos[sp].data.ac, Int.(model.individuals.phytos[sp].data.idx), 
                    false, model.arch)

        ##### copy back to model.individuals.phytos[sp].data
        zero_tmp!(model.individuals.phytos[sp].data)
        get_tind!(model.timestepper.tmp.idx, model.timestepper.tmp.ac)
        copyto_tmp!(model.timestepper.tmp, model.individuals.phytos[sp].data, 
                    model.timestepper.tmp.ac, Int.(model.timestepper.tmp.idx), 
                    false, model.arch)
    end


    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end