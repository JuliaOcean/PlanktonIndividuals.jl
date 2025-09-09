##### update physiological attributes of each individual
function plankton_growth!(plank, trs, rnd, p, ΔT, t, arch::Architecture)
    calc_trs_uptake!(plank, trs, p, ΔT, arch)
    update_state_1!(plank, ΔT, arch)

    calc_PS!(plank, trs, p, arch)
    calc_repiration!(plank,trs,p, ΔT, arch)
    calc_carbon_fixation!(plank, trs, p, arch)
    calc_NO3_reduction!(plank, trs, p, ΔT, arch)
    calc_nitrogen_fixation!(plank, trs, p, arch)
    energy_allocation!(plank, p, arch) 

    update_state_2!(plank, ΔT, arch)
    
    calc_ρChl!(plank, trs.par, p, arch)
    calc_BS!(plank, p, arch)
    calc_iron_fluxes!(plank, trs, p, ΔT, arch)
    update_biomass!(plank, p, ΔT, arch)
    update_cellsize!(plank, p, arch)
    update_tdark!(plank, trs, ΔT, arch)

    ##### probabilities of grazing, mortality, and cell division
    calc_graz_quadratic!(plank, trs, p.grz_P, arch)
    calc_mort!(plank, p, arch)
    calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
    get_probability!(plank, rnd, ΔT, arch)
end