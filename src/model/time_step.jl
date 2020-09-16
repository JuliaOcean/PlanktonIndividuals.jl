"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_Struct, ΔT, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1

    model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    ##### clear operating array for physiological calculations
    zero_fields!(model.timestepper.plk)
    model.timestepper.chl .= 0.0

    # for plank in model.individuals.phytos
    #     plankton_advectionRK4!(plank.data, model.arch, model.grid,
    #                           model.velocities.vel₀, model.velocities.vel½, model.velocities.vel₁, ΔT)

    #     plankton_diffusion!(plank.data, model.arch, model.grid, model.params["κhP"], ΔT)

    #     ##### calculate accumulated chla quantity (not concentration)
    #     # acc_chla_field!(model.timestepper.chl, plank.data, Int.(plank.data[:,13:15]), model.arch)
    # end

    ##### calculate PAR
    # calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
    #           model.grid, model.params["kc"], model.params["kw"])

    # plankton_update!(model.individuals.phytos, model.timestepper.ope, model.timestepper.plk,
    #                  model.timestepper.par, model.timestepper.chl, model.diags, model.arch,
    #                  model.input.temp[:,:,:,clock], model.input.PAR[:,:,clock],
    #                  model.nutrients.DOC, model.nutrients.NH4, model.nutrients.NO3, model.nutrients.PO4,
    #                  model.grid, model.params, ΔT, model.t)

    zero_fields!(model.timestepper.Gcs)
    zero_fields!(model.timestepper.MD1)
    zero_fields!(model.timestepper.MD2)
    zero_fields!(model.timestepper.MD3)

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
end

