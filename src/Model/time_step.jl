"""
    TimeStep!(model::PlanktonModel, ΔT::Int64, diags::PlanktonDiagnostics, resultpath::String)
Update physiology processes and tracer field of `PlanktonModel` one time step forward.

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

    ##### abiotic particle advection, diffusion, and update
    for sa in keys(model.individuals.abiotics)
        particles_from_bcs!(model.individuals.abiotics[sa].data, model.timestepper.tracer_temp.DFe.data, 
                            model.individuals.abiotics[sa].bc, model.timestepper.rnd_3d, model.individuals.abiotics[sa].p, 
                            ΔT, model.iteration, model.grid, model.t, model.arch)
        ##### RK4
        particle_advection!(model.individuals.abiotics[sa].data, model.timestepper.velos, model.grid, 
                            model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        particle_diffusion!(model.individuals.abiotics[sa].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], model.bgc_params["κhP"], model.bgc_params["κvP"],
                            ΔT, model.grid, model.arch)
        
        ##### Update
        find_inds!(model.individuals.abiotics[sa].data, model.grid, model.arch)
    end # abiotic particles

    ##### phytoplankton advection, diffusion, and physiological update
    if model.bgc_params["shared_graz"] == 1.0f0 # shared grazing
        @inbounds model.timestepper.pop .= 0.0f0
        for sp in keys(model.individuals.phytos)
            ##### RK4
            particle_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                                model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
            ##### Diffusion
            particle_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                                model.bgc_params["κhP"], model.bgc_params["κhP"], model.bgc_params["κvP"],
                                ΔT, model.grid, model.arch)

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
            calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, 
                      model.timestepper.PARF, model.grid, model.bgc_params["kc"], 
                      model.bgc_params["kw"], ki)
        end
        for sp in keys(model.individuals.phytos)
            find_NPT!(model.timestepper.trs, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.tracers.NH4.data,
                      model.tracers.NO3.data, model.tracers.PO4.data, model.tracers.DOC.data,
                      model.tracers.DFe.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)
            
            plankton_update!(model.individuals.phytos[sp], model.timestepper.trs,
                             model.timestepper.rnd, model.timestepper.plk, 
                             diags.phytos[sp], ΔT, model.t, model.arch, model.mode)
        end
    else # model.bgc_params["shared_graz"] ≠ 1.0 - species-specific grazing
        for sp in keys(model.individuals.phytos)
            ##### RK4
            particle_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                                model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
            ##### Diffusion
            particle_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                                model.bgc_params["κhP"], model.bgc_params["κhP"], model.bgc_params["κvP"],
                                ΔT, model.grid, model.arch)

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
            @inbounds model.timestepper.pop .= 0.0f0
            acc_counts!(model.timestepper.pop, model.individuals.phytos[sp].data.ac,
                        model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi,
                        model.individuals.phytos[sp].data.zi, model.arch)

            find_NPT!(model.timestepper.trs, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.tracers.NH4.data,
                      model.tracers.NO3.data, model.tracers.PO4.data, model.tracers.DOC.data,
                      model.tracers.DFe.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)

            plankton_update!(model.individuals.phytos[sp], model.timestepper.trs,
                                model.timestepper.rnd, model.timestepper.plk, 
                                diags.phytos[sp], ΔT, model.t, model.arch, model.mode)
        end
    end # phytoplankton

    ##### particle-particle interaction
    for pair in model.timestepper.palat.intac
        plank = model.individuals.phytos[pair[1]].data
        plank_p = model.individuals.phytos[pair[1]].p
        abiotic = model.individuals.abiotics[pair[2]].data
        abio_p = model.individuals.abiotics[pair[2]].p
        particle_interaction!(abiotic, plank, plank_p, model.timestepper.intac, abio_p,
                              model.timestepper.rnd, model.grid, model.max_candidates, model.arch)
    end

    ##### particle-particle release
    for pair in model.timestepper.palat.release
        plank = model.individuals.phytos[pair[1]].data
        abiotic = model.individuals.abiotics[pair[2]].data
        abio_p = model.individuals.abiotics[pair[2]].p
        particle_release!(plank, abiotic, model.timestepper.trs, model.timestepper.rnd,
                          abio_p, ΔT, model.t, model.arch)
    end

    ##### diagnostics of particle-particle interaction
    if isempty(model.individuals.abiotics) == false   
        for sp in keys(model.individuals.phytos)
            diags_proc!(diags.phytos[sp].ptc, 
                        model.individuals.phytos[sp].data.ptc, 
                        model.individuals.phytos[sp].data.ac, 
                        model.individuals.phytos[sp].data.xi, 
                        model.individuals.phytos[sp].data.yi, 
                        model.individuals.phytos[sp].data.zi, model.arch)
        end
    end

    for sa in keys(model.individuals.abiotics)
        diags_proc!(diags.abiotics[sa].num, 
                    model.individuals.abiotics[sa].data.ac, 
                    model.individuals.abiotics[sa].data.ac, 
                    model.individuals.abiotics[sa].data.xi, 
                    model.individuals.abiotics[sa].data.yi, 
                    model.individuals.abiotics[sa].data.zi, model.arch)
    end
    
    ##### tracers update
    tracer_update!(model.tracers, model.timestepper.Gcs, model.timestepper.tracer_temp, 
                   model.timestepper.flux_sink, model.arch,
                   model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT, 
                   model.iteration)

    ##### diagnostics for tracers
    @inbounds diags.tracer.PAR .+= model.timestepper.par
    for key in keys(diags.tracer)
        if key in keys(model.tracers)
            @inbounds diags.tracer[key] .+= model.tracers[key].data
        end
    end # tracers

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
    @inbounds model.timestepper.par₀ .= model.timestepper.par

    return nothing
end
