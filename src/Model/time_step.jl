"""
    TimeStep!(model::PlanktonModel, ΔT::Int64, diags::PlanktonDiagnostics, resultpath::String)
Update physiology processes and nutrient field of `PlanktonModel` one time step forward.

Keyword Arguments
=================
- `model`: `PlanktonModel` to be updated one time step forward.
- `ΔT`: The length of a time step.
- `diags`: `PlanktonDiagnostics` to be updated.
- `resultpath` (optional): The file path to store model output. 
"""
function TimeStep!(model::PlanktonModel, ΔT, diags::PlanktonDiagnostics, resultspath::String)
    model.t = model.t+ΔT
    model.iteration = model.iteration+1

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.Chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing
    ##### plankton advection and diffusion
    for sp in keys(model.individuals.phytos)
        ##### RK4
        plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### AB2 only for 1 species setup
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, 0.1, model.timestepper.vel₁, ΔT, model.arch)
        ##### Euler-forward
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], ΔT, model.grid, model.arch)

        #### calculate accumulated Chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.Chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.Chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, model.timestepper.PARF,
              model.grid, model.bgc_params["kc"], model.bgc_params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.timestepper.temp, model.timestepper.pop, model.arch)

        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                             model.individuals.phytos[sp].proc, model.individuals.phytos[sp].p,
                             model.timestepper.plk, diags.plankton[sp], ΔT, model.t, model.arch, model.mode)
    end
    write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    @inbounds diags.tracer.PAR .+= model.timestepper.par
    for key in keys(diags.tracer)
        if key in keys(model.nutrients)
            @inbounds diags.tracer[key] .+= model.nutrients[key].data
        end
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.nut_temp, model.arch,
                model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT, model.iteration)

    # write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data

    return nothing
end

function TimeStep!(model::PlanktonModel, ΔT, ::Nothing, ::Nothing)
    model.t = model.t+ΔT
    model.iteration = model.iteration+1

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.Chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing
    ##### plankton advection and diffusion
    for sp in keys(model.individuals.phytos)
        ##### RK4
        plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], ΔT, model.grid, model.arch)

        #### calculate accumulated Chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.Chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.Chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, model.timestepper.PARF,
              model.grid, model.bgc_params["kc"], model.bgc_params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.timestepper.temp, model.timestepper.pop, model.arch)

        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                             model.individuals.phytos[sp].proc, model.individuals.phytos[sp].p,
                             model.timestepper.plk, nothing, ΔT, model.t, model.arch, model.mode)
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.nut_temp, model.arch,
                model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT, model.iteration)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data

    return nothing
end