"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_Struct, ΔT, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1
    diag_t = model.t ÷ model.params["diag_freq"] + 1

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0

    ##### plankton advection
    for plank in model.individuals.phytos
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(plank.data, model.timestepper.rnd, model.params["κhP"], ΔT)
        periodic_domain!(plank.data, plank.data.ac, model.grid)

        plankton_advectionRK4!(plank.data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT)

        # plankton_advection!(plank.data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT)

        ##### calculate accumulated chla quantity (not concentration)
        find_inds!(plank.data, plank.data.ac, model.grid)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    plank.data.chl, plank.data.ac, Int.(plank.data.xi),
                    Int.(plank.data.yi), Int.(plank.data.zi), model.grid, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for plank in model.individuals.phytos
        gen_rand_plk!(model.timestepper.rnd, model.arch)
        plankton_update!(plank.data, model.timestepper.nuts, plank.proc, model.timestepper.rnd, 
                         model.timestepper.par, model.timestepper.pop, model.input.temp[:,:,:,clock], 
                         model.nutrients, model.grid, plank.p, ΔT, model.t)

        calc_consume!(model.timestepper.plk.DIC.data, model.timestepper.plk.DOC.data, 
                      model.timestepper.plk.NH4.data, model.timestepper.plk.NO3.data, 
                      model.timestepper.plk.PO4.data, plank.proc, plank.data.ac, 
                      Int.(plank.data.xi), Int.(plank.data.yi), Int.(plank.data.zi), ΔT, model.grid, model.arch)
        # ##### diagnostics for each species and grazing
        # diags!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        ##### grazing
        zero_tmp!(model.timestepper.tmp)
        grazing!(plank.data, model.timestepper.tmp, model.arch, model.grid, model.timestepper.plk, plank.p)

        ###### mortality and its diagnostic
        # diags_mort!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        zero_tmp!(model.timestepper.tmp)
        mortality!(plank.data, model.timestepper.tmp, model.arch, model.grid, model.timestepper.plk, plank.p)

        ###### cell division and its diagnostic
        # diags_dvid!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        ##### division
        zero_tmp!(model.timestepper.tmp)
        dvidnum = dot(plank.data.dvid, plank.data.ac)
        divide!(plank.data, model.timestepper.tmp, dvidnum, model.arch)

        ##### tidy up plank.data, copy to tmp, to the end of divided individuals
        zero_tmp!(model.timestepper.tmp)
        get_tind!(plank.data.idx, plank.data.ac)
        plank.data.idx .= plank.data.idx .+ dvidnum*2
        copyto_tmp!(plank.data, model.timestepper.tmp, plank.data.ac, Int.(plank.data.idx), false, model.arch)

        ##### copy back to plank.data
        zero_tmp!(plank.data)
        get_tind!(model.timestepper.tmp.idx, model.timestepper.tmp.ac)
        copyto_tmp!(model.timestepper.tmp, plank.data, model.timestepper.tmp.ac, 
                    Int.(model.timestepper.tmp.idx), false, model.arch)

    end
    write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    # @inbounds model.diags.tr[:,:,:,diag_t,1] .+= model.timestepper.par
    # @inbounds model.diags.tr[:,:,:,diag_t,2] .+= interior(model.nutrients.NO3.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,3] .+= interior(model.nutrients.NH4.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,4] .+= interior(model.nutrients.PO4.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,5] .+= interior(model.nutrients.DOC.data, model.grid)

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end
