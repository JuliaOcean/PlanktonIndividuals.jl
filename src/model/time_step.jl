"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology processes and nutrient field of `PI_Model` one time step forward.

Keyword Arguments
=================
- `model`: `PI_Model` to be updated one time step forward
- `ΔT`: The length of a time step
- `resultpath` (optional): The file path to store model output. 
"""
function PI_TimeStep!(model::PI_Model, ΔT, resultspath::String)
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
        gen_rand!(model.timestepper.rnd, model.arch)
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd, model.params["κhP"], ΔT, model.arch)
        periodic_domain!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid, model.arch)

        plankton_advectionRK4!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)

        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)

        #### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.input.temp[:,:,:,clock], model.timestepper.pop,
                  model.individuals.phytos[sp].p)

        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts, 
                         model.individuals.phytos[sp].proc, model.timestepper.rnd, 
                         model.individuals.phytos[sp].p, ΔT, model.t, model.arch)

        calc_consume!(model.timestepper.plk.DIC.data, model.timestepper.plk.DOC.data, 
                      model.timestepper.plk.NH4.data, model.timestepper.plk.NO3.data, 
                      model.timestepper.plk.PO4.data, model.individuals.phytos[sp].proc, 
                      model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      ΔT, model.arch)
        ##### diagnostics of processes for each species
        diags_spcs!(model.diags.spcs[sp], model.individuals.phytos[sp].proc, model.individuals.phytos[sp].data,
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, model.arch)
        ##### check the probabilities every 10 mins
        if model.t%300 == 1
            ##### grazing and its diagnostic
            diags_proc!(model.diags.spcs[sp].graz, model.individuals.phytos[sp].data.graz,
                        model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                        model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                        model.arch)

            grazing!(model.individuals.phytos[sp].data, model.arch, 
                     model.timestepper.plk, model.individuals.phytos[sp].p)

            ###### mortality and its diagnostic
            diags_proc!(model.diags.spcs[sp].mort, model.individuals.phytos[sp].data.mort,
                        model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                        model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                        model.arch)

            mortality!(model.individuals.phytos[sp].data, model.arch, 
                       model.timestepper.plk, model.individuals.phytos[sp].p)

            ###### cell division diagnostic
            diags_proc!(model.diags.spcs[sp].dvid, model.individuals.phytos[sp].data.dvid,
                        model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                        model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                        model.arch)

            ##### division
            ##### check if the number of individuals exceeded
            dvidnum = dot(model.individuals.phytos[sp].data.dvid, model.individuals.phytos[sp].data.ac)
            deactive_ind = findall(x -> x == 0.0, model.individuals.phytos[sp].data.ac)
            if dvidnum > length(deactive_ind)
                throw(ArgumentError("number of individual exceeds the capacity $(model.params["λ"])N"))
            end
            ##### do not copy inactive individuals
            model.individuals.phytos[sp].data.dvid .*= model.individuals.phytos[sp].data.ac
            divide!(model.individuals.phytos[sp].data, deactive_ind, model.arch)
        end

        ##### diagnostic for individual distribution
        diags_proc!(model.diags.spcs[sp].num, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                    model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                    model.arch)

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

    return nothing
end

function PI_TimeStep!(model::PI_Model, ΔT)
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
        gen_rand!(model.timestepper.rnd, model.arch)
        gen_rand_adv!(model.timestepper.rnd, model.arch)
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd, model.params["κhP"], ΔT, model.arch)
        periodic_domain!(model.individuals.phytos[sp].data, model.individuals.phytos[sp].data.ac, model.grid, model.arch)

        plankton_advectionRK4!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)

        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)

        #### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.input.PARF[:,:,clock],
              model.grid, model.params["kc"], model.params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.input.temp[:,:,:,clock], model.timestepper.pop,
                  model.individuals.phytos[sp].p)

        plankton_update!(model.individuals.phytos[sp].data, model.timestepper.nuts, 
                         model.individuals.phytos[sp].proc, model.timestepper.rnd, 
                         model.individuals.phytos[sp].p, ΔT, model.t, model.arch)

        calc_consume!(model.timestepper.plk.DIC.data, model.timestepper.plk.DOC.data, 
                      model.timestepper.plk.NH4.data, model.timestepper.plk.NO3.data, 
                      model.timestepper.plk.PO4.data, model.individuals.phytos[sp].proc, 
                      model.individuals.phytos[sp].data.ac, model.individuals.phytos[sp].data.xi, 
                      model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi,
                      ΔT, model.arch)
        ##### check the probabilities every 10 mins
        if model.t%300 == 1
            ##### grazing
            grazing!(model.individuals.phytos[sp].data, model.arch, 
                     model.timestepper.plk, model.individuals.phytos[sp].p)

            ###### mortality
            mortality!(model.individuals.phytos[sp].data, model.arch, 
                       model.timestepper.plk, model.individuals.phytos[sp].p)

            ##### division
            ##### check if the number of individuals exceeded
            dvidnum = dot(model.individuals.phytos[sp].data.dvid, model.individuals.phytos[sp].data.ac)
            deactive_ind = findall(x -> x == 0.0, model.individuals.phytos[sp].data.ac)
            if dvidnum > length(deactive_ind)
                throw(ArgumentError("number of individual exceeds the capacity $(model.params["λ"])N"))
            end
            ##### do not copy inactive individuals
            model.individuals.phytos[sp].data.dvid .*= model.individuals.phytos[sp].data.ac
            divide!(model.individuals.phytos[sp].data, deactive_ind, model.arch)
        end
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.MD1,
                model.timestepper.MD2, model.timestepper.MD3, model.arch,
                model.grid, model.params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
    return nothing
end