##### update physiological attributes of each individual
function plankton_growth!(plank, trs, rnd, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, trs, p, ΔT, arch)

    update_quotas_1!(plank, ΔT, arch)

    calc_organic_uptake!(plank, trs, p, ΔT, arch)

    calc_ρChl!(plank, trs.par, p, arch)

    calc_respir!(plank, trs.T, p, arch)

    update_quotas_2!(plank, ΔT, p, arch)

    calc_BS!(plank, p, arch)
    update_biomass!(plank, p, ΔT, arch)
    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    calc_graz_quadratic!(plank, trs, p.grz_P, arch)
    calc_mort!(plank, p, arch)
    calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
    get_probability!(plank, rnd, ΔT, arch)
end