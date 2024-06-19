##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch::Architecture)
    calc_PS!(plank, nuts, p, arch)
    update_state_1!(plank, ΔT, arch)
    calc_carbon_fluxes!(plank, nuts, p, arch)
    calc_nuts_uptake!(plank, nuts, p, arch)
    update_state_2!(plank, ΔT, arch)
    calc_NO3_reduc!(plank, nuts, p, arch)
    update_state_3!(plank, ΔT, arch)
    calc_ρChl!(plank, nuts.par, p, arch)
    calc_BS!(plank, p, arch)
    update_biomass!(plank, p, ΔT, arch)
    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10.0f0,3600.0f0÷ΔT))) == 0.0f0 
        calc_graz_quadratic!(plank, nuts, p.grz_P, arch)
        calc_mort!(plank, p, arch)
        calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, rnd, ΔT, arch)

        ##### thermal mortality (WIP)
        if p.ther_mort == 1.0f0
            # thermal_mort!(plank, nuts, p)
        end
    else
        @inbounds plank.graz .= 0.0f0
        @inbounds plank.mort .= 0.0f0
        @inbounds plank.dvid .= 0.0f0
    end

end