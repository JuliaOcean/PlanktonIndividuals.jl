##### update physiological attributes of each individual
function plankton_growth!(plank, trs, rnd, p, ΔT, t, g::AbstractGrid, arch::Architecture)

    calc_PS!(plank, trs, p, arch)

    calc_respiration!(plank, trs.T, p, ΔT, arch)

    calc_inorganic_uptake!(plank, trs, p, ΔT, g, arch)

    calc_uptake_energy_alloc!(plank, p, arch)

    update_quotas_1!(plank, ΔT, arch)

    calc_carbon_fixation!(plank, trs.T, p, arch)

    calc_NO3_reduction!(plank, trs.T, p, ΔT, arch)

    calc_nitrogen_fixation!(plank, trs.T, p, arch::Architecture)

    calc_CN_energy_alloc!(plank, p, arch)

    calc_organic_uptake!(plank, trs, p, ΔT, arch)

    update_quotas_2!(plank, ΔT, p, arch)

    calc_Chl_synthesis!(plank, trs.par, p, arch)

    calc_BS!(plank, trs, p, arch, ΔT)
    
    calc_BS_energy_alloc!(plank, p, arch)

    update_biomass!(plank, p, ΔT, arch)

    calc_exudation!(plank, p, arch)

    update_CH!(plank, arch)

    calc_degradation!(plank, p, trs, ΔT, arch)

    update_quotas_3!(plank, ΔT, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    calc_graz_quadratic!(plank, trs, p.grz_P, arch)
    calc_MM_mort!(plank, p, arch)
    ##### Bernouli-like distribution
    calc_MM_dvid!(plank, p, arch)
    get_probability!(plank, rnd, ΔT, arch)
end