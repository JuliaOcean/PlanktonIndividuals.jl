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
function TimeStep!(model::PlanktonModel, ΔT, diags::PlanktonDiagnostics)
    # model.t = model.t+ΔT
    model.iteration = model.iteration+1
    model.t = model.iteration * ΔT 

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5f0
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5f0
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5f0

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.Chl .= 0.0f0
    @inbounds model.timestepper.pop .= 0.0f0

    ##### abiotic particle advection, diffusion, and update
    for sp in keys(model.individuals.abiotics)
        ##### RK4
        particle_advection!(model.individuals.abiotics[sp].data, model.timestepper.velos, model.grid, 
                            model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        particle_diffusion!(model.individuals.abiotics[sp].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], ΔT, model.grid, model.arch)
        
        ##### Update
        find_inds!(model.individuals.abiotics[sp].data, model.grid, model.arch)
        particle_update!(model.individuals.abiotics[sp].data, model.nutrients.CHO.data,
                         model.timestepper.plk.CHO.data, model.individuals.abiotics[sp].p,
                         diags.abiotics[sp], ΔT, model.acch)
        
    end # abiotic particles

    ##### phytoplankton advection, diffusion, and physiological update
    if model.bgc_params["shared_graz"] == 1.0f0 # shared grazing
        for sp in keys(model.individuals.phytos)
            ##### RK4
            particle_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                                model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
            ##### Diffusion
            particle_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                                model.bgc_params["κhP"], ΔT, model.grid, model.arch)

            #### calculate accumulated Chla quantity (not concentration) and population
            find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
            acc_chl!(model.timestepper.Chl, model.individuals.phytos[sp].data.Chl,
                     model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi,
                     model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, model.arch)
            acc_counts!(model.timestepper.pop, model.individuals.phytos[sp].data.ac,
                        model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi,
                        model.individuals.phytos[sp].data.zi, model.arch)
        end
        ##### calculate PAR
        for ki in 1:model.grid.Nz
            calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, model.timestepper.PARF,
                      model.grid, model.bgc_params["kc"], model.bgc_params["kw"], ki)
        end
        for sp in keys(model.individuals.phytos)
            find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data,
                      model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                      model.nutrients.FeT.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)
            
            plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                             model.timestepper.rnd, model.individuals.phytos[sp].p,
                             model.timestepper.plk, diags.plankton[sp], ΔT, model.t, model.arch, model.mode)
        end
    else # model.bgc_params["shared_graz"] ≠ 1.0 - species-specific grazing
        for sp in keys(model.individuals.phytos)
            ##### RK4
            particle_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                                model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
            ##### Diffusion
            particle_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                                model.bgc_params["κhP"], ΔT, model.grid, model.arch)

            #### calculate accumulated Chla quantity (not concentration) and population
            find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
            acc_chl!(model.timestepper.Chl, model.individuals.phytos[sp].data.Chl,
                     model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi,
                     model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, model.arch)
        end
        ##### calculate PAR
        for ki in 1:model.grid.Nz
            calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, model.timestepper.PARF,
                      model.grid, model.bgc_params["kc"], model.bgc_params["kw"], ki)
        end
        for sp in keys(model.individuals.phytos)
            acc_counts!(model.timestepper.pop, model.individuals.phytos[sp].data.ac,
                        model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi,
                        model.individuals.phytos[sp].data.zi, model.arch)

            find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data,
                      model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                      model.nutrients.FeT.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)

            plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                                model.timestepper.rnd, model.individuals.phytos[sp].p,
                                model.timestepper.plk, diags.plankton[sp], ΔT, model.t, model.arch, model.mode)
            @inbounds model.timestepper.pop .= 0.0f0
        end
    end # phytoplankton

    ##### diagnostics for nutrients
    @inbounds diags.tracer.PAR .+= model.timestepper.par
    @inbounds diags.tracer.T .+= model.timestepper.temp
    for key in keys(diags.tracer)
        if key in keys(model.nutrients)
            @inbounds diags.tracer[key] .+= model.nutrients[key].data
        end
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.nut_temp, model.arch,
                model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT, model.iteration)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
    @inbounds model.timestepper.par₀ .= model.timestepper.par

    return nothing
end
