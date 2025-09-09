##### update physiological attributes of each individual
function plankton_growth!(plank, trs, rnd, p, ΔT, t, arch::Architecture)
    calc_PS!(plank, trs, p, arch)

    calc_BS!(plank, p, arch)

    calc_repair!(plank, p, arch)

    calc_thermal_damage!(plank, trs.T, p, ΔT, arch)

    calc_respir!(plank, trs.T, p, arch)

    update_quotas!(plank, ΔT, arch)

    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    calc_graz_quadratic!(plank, trs, p.grz_P, arch)
    calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)

    ##### thermal mortality
    if p.thermal == 1.0f0
        calc_thermal_mort!(plank, p, arch)
    else
        calc_mort!(plank, p, arch)
    end

    get_probability!(plank, rnd, ΔT, arch)
end
