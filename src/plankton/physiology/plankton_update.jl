##### update physiological attributes of each individual
function plankton_update!(plank, nuts, proc, rnd, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, proc, nuts, p, arch)

    update_quotas_1!(plank, proc, ΔT, arch)

    calc_organic_uptake!(plank, proc, nuts, p, arch)

    calc_ρchl!(plank, proc, nuts, p.Chl2N, arch)

    calc_respir!(plank, proc, nuts, p, arch)

    update_quotas_2!(plank, proc, ΔT, p, arch)

    calc_BS!(plank, proc, p, arch)
    update_biomass!(plank, proc, p, ΔT, arch)
    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 mins
    if t%600 == 0
        calc_graz_quadratic!(plank, proc, nuts, p.grz_P, arch)
        calc_mort!(plank, proc, p, arch)
        calc_dvid!(plank, proc, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, proc, arch)
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end
end