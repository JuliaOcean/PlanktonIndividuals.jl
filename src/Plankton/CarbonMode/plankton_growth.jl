##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, nuts, p, arch)

    calc_respir!(plank, nuts.T, p, arch)

    update_quotas!(plank, ΔT, arch)

    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
        calc_graz_quadratic!(plank, nuts, p.grz_P, arch)
        calc_mort!(plank, p, arch)
        calc_dvid!(plank, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, rnd, ΔT, arch)

        ##### thermal mortality
        if p.ther_mort == 1
            thermal_mort!(plank, nuts, p)
        end
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end

end