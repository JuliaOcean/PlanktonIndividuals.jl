"""
    PI_TimeStep!(model, RunParam, velᵇ)
Update physiology part and nutrient field of 'model' one time step forward
"""
function PI_TimeStep!(model::Model_Struct, ΔT, vel₁::NamedTuple, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1

    model.velocities.vel₁ = vel₁
    model.velocities.vel½ = (u = (model.velocities.vel₀.u .+ model.velocities.vel₁.u) .* 0.5,
                             v = (model.velocities.vel₀.v .+ model.velocities.vel₁.v) .* 0.5,
                             w = (model.velocities.vel₀.w .+ model.velocities.vel₁.w) .* 0.5)

    for plank in model.individuals.phytos
        plankton_advectionRK4!(plank.data, model.arch, model.grid,
                              model.velocities.vel₀, model.velocities.vel½, model.velocities.vel₁, ΔT)

        plankton_diffusion!(plank.data, model.arch, model.grid, model.params["κhP"], ΔT)
    end

    ##### clear operating array for physiological calculations
    zero_fields!(model.timestepper.plk)

    # plankton_update!(model.individuals.phytos, model.timestepper.ope, model.timestepper.plk,
    #                  model.timestepper.par, model.timestepper.chl, model.diags, model.arch,
    #                  model.input.temp[:,:,:,clock], model.input.PAR[:,:,clock],
    #                  model.nutrients.DOC, model.nutrients.NH4, model.nutrients.NO3, model.nutrients.PO4,
    #                  model.grid, model.params, ΔT, model.t)


    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    model.velocities.vel₀ = vel₁
end
function PI_TimeStep!(model::Model_Struct, ΔT, vel₁::NamedTuple)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1
    consume = nutrients_init(model.arch, model.grid)
    model.velocities.vel₁ = vel₁

    plankton_advectionRK4!(model.individuals.phytos, model.arch, model.grid,
                           model.velocities.vel₀, model.velocities.vel₁, ΔT)
    plankton_diffusion!(model.individuals.phytos, model.arch, model.grid,
                        model.params["κhP"], ΔT)

    plankton_update!(model.individuals.phytos, consume, model.diags, model.arch,
                     model.input.temp[:,:,:,clock], model.input.PAR[:,:,clock],
                     model.nutrients.DOC, model.nutrients.NH4, model.nutrients.NO3, model.nutrients.PO4,
                     model.grid, model.params, ΔT, model.t)

    gtr = nutrients_init(model.arch, model.grid)

    nut_update!(model.nutrients, gtr, model.arch, model.grid, model.params, vel₁, consume, ΔT)

    model.velocities.vel₀ = vel₁
end
function PI_TimeStep!(model::Model_Struct, ΔT, resultspath::String)
    model.t = model.t+ΔT
    clock = model.t % 86400 ÷ ΔT + 1
    consume = nutrients_init(model.arch, model.grid)

    plankton_update!(model.individuals.phytos, consume, model.diags, model.arch,
                     model.input.temp[:,:,:,clock], model.input.PAR[:,:,clock],
                     model.nutrients.DOC, model.nutrients.NH4, model.nutrients.NO3, model.nutrients.PO4,
                     model.grid, model.params, ΔT, model.t)

    gtr = nutrients_init(model.arch, model.grid)

    nut_update!(model.nutrients, gtr, model.arch, model.grid, model.params, consume, ΔT)

    write_nut_cons(model.grid, gtr, model.nutrients, model.t, resultspath)
end

