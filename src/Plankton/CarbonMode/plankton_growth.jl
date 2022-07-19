##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, nuts, p, arch)

    calc_respir!(plank, nuts.T, p, arch)

    update_quotas!(plank, ΔT, arch)

    update_cellsize!(plank, p, arch)

    calc_thermal_history!(plank, nuts, p, ΔT, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
        calc_graz_quadratic!(plank, nuts, p.grz_P, arch)
        calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)

        ##### thermal mortality (WIP)
        if p.ther_mort == 1
            calc_thermal_mort!(plank, p, arch)
        else
            calc_mort!(plank, p, arch)
        end

        get_probability!(plank, rnd, ΔT, arch)
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end

    # ##### reset thermal history factor every day
    # if t%86400 == 0
    #     plank.Th .= 0.0
    # end

end