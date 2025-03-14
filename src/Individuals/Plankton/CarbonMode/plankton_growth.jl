##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch::Architecture)
    calc_PS!(plank, nuts, p, arch)

    calc_BS!(plank, p, arch)

    calc_repair!(plank, p, arch)

    calc_thermal_damage!(plank, nuts.T, p, ΔT, arch)

    calc_respir!(plank, nuts.T, p, arch)

    update_quotas!(plank, ΔT, arch)

    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10.0f0,3600.0f0÷ΔT))) == 0.0f0
        calc_graz_quadratic!(plank, nuts, p.grz_P, arch)
        calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)

        ##### thermal mortality
        if p.thermal == 1.0f0
            calc_thermal_mort!(plank, p, arch)
        else
            calc_mort!(plank, p, arch)
        end

        get_probability!(plank, rnd, ΔT, arch)
    else
        @inbounds plank.graz .= 0.0f0
        @inbounds plank.mort .= 0.0f0
        @inbounds plank.dvid .= 0.0f0
    end
end
