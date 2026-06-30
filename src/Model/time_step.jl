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

    ##### abiotic particle motion
    for sa in eachindex(model.individuals.abiotics)
        particles_from_bcs!(model.individuals.abiotics[sa].data, model.timestepper.tracer_temp.DFe.data, 
                            model.individuals.abiotics[sa].bc, model.timestepper.rnd_3d, model.individuals.abiotics[sa].p, 
                            ΔT, model.iteration, model.grid, model.t, model.arch)
        ##### particle motion
        particle_motion!(model.individuals.abiotics[sa].data, model.timestepper.velos, model.grid, 
                         model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, 
                         model.timestepper.rnd, model.bgc_params["κhP"], model.bgc_params["κhP"], 
                         model.bgc_params["κvP"],ΔT, model.arch)
    end # abiotic particles

    ##### phytoplankton motion
    for sp in eachindex(model.individuals.phytos)
        particle_motion!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid, 
                         model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, 
                         model.timestepper.rnd, model.bgc_params["κhP"], model.bgc_params["κhP"], 
                         model.bgc_params["κvP"],ΔT, model.arch)
    end # phytoplankton motion

    ##### colonies motion
    for cl in eachindex(model.individuals.colonies)
        colony_motion!(model.individuals.colonies[cl].spcs, model.timestepper.velos, model.grid, 
                       model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, 
                       model.timestepper.rnd, model.bgc_params["κhP"], model.bgc_params["κhP"], 
                       model.bgc_params["κvP"],ΔT, model.arch)
    end # colony motion

    ##### calculate accumulated Chla quantity (not concentration)
    for sp in eachindex(model.individuals.phytos)
        acc_chl!(model.timestepper.Chl, model.individuals.phytos[sp].data.Chl,
                 model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi,
                 model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, model.arch)
    end
    for cl in eachindex(model.individuals.colonies)
        for sp in eachindex(model.individuals.colonies[cl].spcs)
            acc_chl!(model.timestepper.Chl, model.individuals.colonies[cl].spcs[sp].data.Chl,
                     model.individuals.colonies[cl].spcs[sp].data.ac, 
                     model.individuals.colonies[cl].spcs[sp].data.xi,
                     model.individuals.colonies[cl].spcs[sp].data.yi, 
                     model.individuals.colonies[cl].spcs[sp].data.zi, model.arch)
        end
    end # Chla

    ##### calculate PAR
    for ki in 1:model.grid.Nz
        calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, 
                  model.timestepper.PARF, model.grid, model.bgc_params["kc"], 
                  model.bgc_params["kw"], ki)
    end # PAR

    ##### phytoplankton physiological update
    if model.bgc_params["shared_graz"] == 1.0f0 # shared grazing
        @inbounds model.timestepper.pop .= 0.0f0
        for sp in eachindex(model.individuals.phytos)
            #### calculate population
            acc_counts!(model.timestepper.pop, model.individuals.phytos[sp].data.ac,
                        model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi,
                        model.individuals.phytos[sp].data.zi, model.arch)
        end
        for sp in eachindex(model.individuals.phytos)
            find_NPT!(model.timestepper.trs, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.tracers.NH4.data,
                      model.tracers.NO3.data, model.tracers.PO4.data, model.tracers.DOC.data,
                      model.tracers.DFe.data, model.tracers.O2.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)
            
            plankton_update!(model.individuals.phytos[sp], model.timestepper.trs,
                             model.timestepper.rnd, model.timestepper.plk, 
                             diags.phytos[sp], ΔT, model.t, model.arch, model.mode)
        end
    else # model.bgc_params["shared_graz"] ≠ 1.0 - species-specific grazing
        for sp in eachindex(model.individuals.phytos)
            @inbounds model.timestepper.pop .= 0.0f0
            acc_counts!(model.timestepper.pop, model.individuals.phytos[sp].data.ac,
                        model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi,
                        model.individuals.phytos[sp].data.zi, model.arch)

            find_NPT!(model.timestepper.trs, model.individuals.phytos[sp].data.xi,
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      model.individuals.phytos[sp].data.ac, model.tracers.NH4.data,
                      model.tracers.NO3.data, model.tracers.PO4.data, model.tracers.DOC.data,
                      model.tracers.DFe.data, model.tracers.O2.data, model.timestepper.par, model.timestepper.par₀, 
                      model.timestepper.temp, model.timestepper.pop, model.arch)

            plankton_update!(model.individuals.phytos[sp], model.timestepper.trs,
                                model.timestepper.rnd, model.timestepper.plk, 
                                diags.phytos[sp], ΔT, model.t, model.arch, model.mode)
        end
    end # phytoplankton

    ##### colony physiology update
    for cl in eachindex(model.individuals.colonies)
        @inbounds model.timestepper.pop .= 0.0f0
        acc_counts!(model.timestepper.pop, model.individuals.colonies[cl].spcs.sp1.data.ac,
                    model.individuals.colonies[cl].spcs.sp1.data.xi, 
                    model.individuals.colonies[cl].spcs.sp1.data.yi,
                    model.individuals.colonies[cl].spcs.sp1.data.zi, model.arch)

        find_NPT!(model.timestepper.trs, model.individuals.colonies[cl].spcs.sp1.data.xi,
                  model.individuals.colonies[cl].spcs.sp1.data.yi, 
                  model.individuals.colonies[cl].spcs.sp1.data.zi,
                  model.individuals.colonies[cl].spcs.sp1.data.ac, model.tracers.NH4.data,
                  model.tracers.NO3.data, model.tracers.PO4.data, model.tracers.DOC.data,
                  model.tracers.DFe.data, model.tracers.O2.data, model.timestepper.par, model.timestepper.par₀, 
                  model.timestepper.temp, model.timestepper.pop, model.arch)

        colony_update!(model.individuals.colonies[cl], model.timestepper.trs,
                       model.timestepper.rnd, model.timestepper.plk, 
                       diags.colonies[cl], ΔT, model.t, model.arch, model.mode)
    end

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
        for sp in eachindex(model.individuals.phytos)
            diags_proc!(diags.phytos[sp].ptc, 
                        model.individuals.phytos[sp].data.ptc, 
                        model.individuals.phytos[sp].data.ac, 
                        model.individuals.phytos[sp].data.xi, 
                        model.individuals.phytos[sp].data.yi, 
                        model.individuals.phytos[sp].data.zi, model.arch)
        end
    end

    ##### diagnostics for abiotic particles
    for sa in eachindex(model.individuals.abiotics)
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
    for key in eachindex(diags.tracer)
        if key in eachindex(model.tracers)
            @inbounds diags.tracer[key] .+= model.tracers[key].data
        end
    end # tracers

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
    @inbounds model.timestepper.par₀ .= model.timestepper.par

    return nothing
end
