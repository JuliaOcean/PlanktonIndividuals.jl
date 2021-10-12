##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, proc, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, proc, nuts, p, arch)

    calc_respir!(plank, proc, nuts.T, p, arch)

    update_quotas!(plank, proc, ΔT, arch)

    update_cellsize!(plank, p, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
        calc_graz_quadratic!(plank, proc, nuts, p.grz_P, arch)
        calc_mort!(plank, proc, p, arch)
        calc_dvid!(plank, proc, divide_type(p.dvid_type), p, t, arch)
        get_probability!(plank, proc, ΔT, arch)
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end
end