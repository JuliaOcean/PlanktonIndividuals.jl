##### update physiological attributes of each individual
function colony_plankton_growth!(plank, trs, p, ΔT, arch::Architecture)
    calc_NPFe_uptake!(plank, trs, p, ΔT, arch)
    calc_PS!(plank, trs, p, arch)
    calc_carbon_fixation!(plank, trs, p, arch)
    calc_NO3_reduction!(plank, trs, p, ΔT, arch)
    calc_nitrogen_fixation!(plank, trs, p, arch)
    calc_O2_production!(plank, p, arch)
    calc_O2_diffusion!(plank, trs, p, ΔT, arch)
    calc_repiration!(plank,trs, p, ΔT, arch)
    energy_redox_allocation!(plank, p, ΔT, arch)

    calc_ρChl!(plank, trs.par, p, arch)
    calc_BS!(plank, p, arch)
    calc_iron_fluxes!(plank, trs, p, ΔT, arch)
    update_states!(plank, p, ΔT, arch)
    update_cellsize!(plank, p, arch)
    update_tdark!(plank, trs, ΔT, arch)
end

function calc_plankton_population_dynamics(plank, trs, rnd, p, ΔT, t, arch::Architecture)
    ##### probabilities of grazing, mortality, and cell division
    calc_graz_quadratic!(plank, trs, p.grz_P, arch)
    calc_mort!(plank, p, arch)
    calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
    get_probability!(plank, rnd, ΔT, arch)
end