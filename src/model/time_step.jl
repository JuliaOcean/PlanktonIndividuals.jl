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

    ##### plankton advection
    for plank in model.individuals.phytos
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(plank.data, model.timestepper.rnd, model.params["κhP"], ΔT)
        periodic_domain!(plank.data, plank.data.ac, model.grid)

        plankton_advectionRK4!(plank.data, model.timestepper.coord, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT)

        # plankton_advection!(plank.data, model.timestepper.coord, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT)

        ##### calculate accumulated chla quantity (not concentration)
        find_inds!(plank.data, model.timestepper.coord, plank.data.ac, model.grid)
        acc_counts!(model.timestepper.cts.chl, model.timestepper.cts.pop,
                    plank.data.chl, plank.data.ac, Int.(model.timestepper.coord.x),
                    Int.(model.timestepper.coord.y), Int.(model.timestepper.coord.z), model.arch)
    end

    ##### calculate PAR
    @inbounds sum!(model.timestepper.chl, model.timestepper.cts.chl)
    @inbounds sum!(model.timestepper.pop, model.timestepper.cts.pop)
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for plank in model.individuals.phytos
        model.timestepper.tmp .= 0.0
        gen_rand_plk!(model.timestepper.rnd, model.arch)
        plankton_update!(plank.data, model.timestepper.nuts, model.timestepper.proc, model.timestepper.coord, 
                         model.timestepper.rnd, model.timestepper.par, model.timestepper.pop, 
                         model.input.temp[:,:,:,clock], model.nutrients, model.grid, plank.p, ΔT, model.t)

        # calc_consume!(model.timestepper.cts.DIC, model.timestepper.cts.DOC, model.timestepper.cts.NH4, 
        #               model.timestepper.cts.NO3, model.timestepper.cts.PO4, model.timestepper.proc, plank.data.ac, 
        #               Int.(model.timestepper.coord.x), Int.(model.timestepper.coord.y), 
        #               Int.(model.timestepper.coord.z), ΔT, model.arch)
        # ##### diagnostics for each species and grazing
        # diags!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        # ##### grazing
        # model.timestepper.tmp .= 0.0
        # grazing!(plank.data, model.timestepper.tmp, model.arch,
        #         model.grid, model.timestepper.plk, plank.p)

        # ###### mortality and its diagnostic
        # diags_mort!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        # model.timestepper.tmp .= 0.0
        # mortality!(plank.data, model.timestepper.tmp, model.arch,
        #           model.grid, model.timestepper.plk, plank.p)

        # ###### cell division and its diagnostic
        # diags_dvid!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, diag_t)

        # ##### tidy up plank.data
        # model.timestepper.tmp .= 0.0
        # CUDA.@allowscalar plank.active_num = floor(Int, sum(plank.data[:,58], dims=1)[1])
        # copyto_tmp!(plank.data, model.timestepper.tmp, plank.data[:,58], Int.(plank.data[:,59]), false, model.arch)
        # plank.data .= copy(model.timestepper.tmp)

        # ##### copy individuals which are ready to divide to the end of plank.data
        # divide!(plank.data, model.arch, plank.num)
        # CUDA.@allowscalar plank.num = floor(Int, sum(plank.data[:,58], dims=1)[1])
    end
    # write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    # @inbounds model.diags.tr[:,:,:,diag_t,1] .+= model.timestepper.par
    # @inbounds model.diags.tr[:,:,:,diag_t,2] .+= interior(model.nutrients.NO3.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,3] .+= interior(model.nutrients.NH4.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,4] .+= interior(model.nutrients.PO4.data, model.grid)
    # @inbounds model.diags.tr[:,:,:,diag_t,5] .+= interior(model.nutrients.DOC.data, model.grid)

    ##### sum up nutrient consumption counts into nutrient tendencies
    # pcm_to_Gcs!(model.timestepper.plk, model.timestepper.cts)

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end

function zero_timestepper_flds!(ts::timestepper)
    zero_fields!(ts.plk)
    @inbounds ts.chl .= 0.0
    @inbounds ts.pop .= 0.0
    @inbounds ts.cts.chl .= 0.0
    @inbounds ts.cts.pop .= 0.0
    # @inbounds ts.cts.DIC .= 0.0
    # @inbounds ts.cts.NH4 .= 0.0
    # @inbounds ts.cts.NO3 .= 0.0
    # @inbounds ts.cts.PO4 .= 0.0
end
