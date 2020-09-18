"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_Struct, ΔT, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1
    diag_t = model.t ÷ model.params["diag_freq"] + 1

    model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    ##### clear operating array for physiological calculations
    zero_fields!(model.timestepper.plk)
    model.timestepper.chl .= 0.0

    for plank in model.individuals.phytos
        plankton_advectionRK4!(plank.data, model.arch, model.grid,
                              model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT)

        gen_rand_adv!(plank.data, model.arch)
        plankton_diffusion!(plank.data, model.arch, model.params["κhP"], ΔT)
        in_domain!(plank.data, model.arch, model.grid)

        ##### calculate accumulated chla quantity (not concentration)
        find_inds!(plank.data, model.arch, model.grid, 12, 0)
        acc_chla_field!(model.timestepper.chl, plank.data, Int.(plank.data[:,13:15]), model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    for plank in model.individuals.phytos
        plank.data[:,58:60] .= rand!(rng_type(model.arch), plank.data[:,58:60])
        plankton_update!(plank.data, model.timestepper.plk,
                         model.timestepper.par, model.arch, model.input.temp[:,:,:,clock],
                         model.nutrients.DOC.data, model.nutrients.NH4.data, model.nutrients.NO3.data,
                         model.nutrients.PO4.data, model.grid, plank.p, ΔT, model.t)

        ###### diagnostics for each species and grazing
        sum_diags!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch,
                   model.grid, diag_t)

        ##### grazing
        plank.data = grazing!(plank.data, model.timestepper.plk, model.arch, model.grid, plank.p)

        ###### mortality and its diagnostic
        sum_diags_mort!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, model.grid, diag_t)
        plank.data = mortality!(plank.data, model.timestepper.plk, model.arch, model.grid, plank.p)

        ###### cell division and its diagnostic
        sum_diags_dvid!(model.diags.spcs, plank.data, Int.(plank.data[:,13:15]), plank.sp, model.arch, model.grid, diag_t)
        plank.data = divide!(plank.data)
    end
    write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    model.diags.tr[:,:,:,diag_t,1] += model.timestepper.par
    model.diags.tr[:,:,:,diag_t,2] += interior(model.nutrients.NO3.data, model.grid)
    model.diags.tr[:,:,:,diag_t,3] += interior(model.nutrients.NH4.data, model.grid)
    model.diags.tr[:,:,:,diag_t,4] += interior(model.nutrients.PO4.data, model.grid)
    model.diags.tr[:,:,:,diag_t,5] += interior(model.nutrients.DOC.data, model.grid)

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end

