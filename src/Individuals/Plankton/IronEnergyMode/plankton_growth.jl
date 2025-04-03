##### update physiological attributes of each individual
function plankton_growth!(plank, trs, rnd, p, ΔT, t, arch::Architecture)
    calc_PS!(plank, trs, p, arch)
    calc_repiration!(plank,trs,p, ΔT, arch)
    update_state_1!(plank, ΔT, arch)
    calc_carbon_fixation!(plank, trs, p, arch)
    calc_trs_uptake!(plank, trs, p, ΔT, arch)
    update_state_2!(plank, ΔT, arch)
    calc_NO3_reduction!(plank, trs, p, ΔT, arch)
    calc_nitrogen_fixation!(plank, trs, p, arch)
    update_state_3!(plank, ΔT, arch)
    calc_fADP!(plank, p, ΔT, arch)
    update_state_4!(plank, ΔT, arch)
    calc_ρChl!(plank, trs.par, p, arch)
    calc_BS!(plank, p, arch)
    calc_iron_fluxes!(plank, trs, p, ΔT, arch)
    update_biomass!(plank, p, ΔT, arch)
    update_cellsize!(plank, p, arch)
    update_tdark!(plank, trs, ΔT, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10.0f0,3600.0f0÷ΔT))) == 0.0f0 
        calc_graz_quadratic!(plank, trs, p.grz_P, arch)
        calc_mort!(plank, p, arch)
        calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, rnd, ΔT, arch)
    else
        @inbounds plank.graz .= 0.0f0
        @inbounds plank.mort .= 0.0f0
        @inbounds plank.dvid .= 0.0f0
    end

end