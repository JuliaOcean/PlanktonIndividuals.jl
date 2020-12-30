##### update physiological attributes of each individual
function plankton_update!(plank, nuts, proc, rnd, p, ΔT, t, arch::Architecture)
    ##### inorganic nutrient uptake
    calc_inorganic_uptake!(plank, proc, nuts, p, arch)

    ##### update C, N, P quotas
    update_quotas_1!(plank, proc, ΔT, arch)

    ##### Organic carbon uptake if necessary
    calc_organic_uptake!(plank, proc, nuts, p, arch)

    ##### Chla
    calc_ρchl!(plank, proc, nuts, p.Chl2N, arch)

    ##### respiration
    calc_respir!(plank, proc, nuts, p, arch)

    ##### update C, N, P quotas
    update_quotas_2!(plank, proc, ΔT, p, arch)

    ##### Biosynthesis
    calc_BS!(plank, proc, p, arch)
    update_biomass!(plank, proc, p, ΔT, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 5 mins
    if t%300 == 1
        gen_rand!(rnd, arch)
        ##### grazing
        calc_graz_quadratic!(plank, proc, nuts, p.grz_P, arch)

        ##### mortality
        calc_mort!(plank, proc, p, arch)

        ##### cell division
        calc_dvid!(plank, proc, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, proc, rnd)
    else
        plank.graz .= 0.0
        plank.mort .= 0.0
        plank.dvid .= 0.0
    end
end