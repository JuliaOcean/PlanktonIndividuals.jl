##### update physiological attributes of each individual
function plankton_growth!(plank, nuts, rnd, p, ΔT, t, arch::Architecture)
    calc_inorganic_uptake!(plank, nuts, p, arch)

    update_quotas_1!(plank, ΔT, arch)

    calc_organic_uptake!(plank, nuts, p, arch)

    calc_ρChl!(plank, nuts.par, p, arch)

    calc_respir!(plank, nuts.T, p, arch)

    update_quotas_2!(plank, ΔT, p, arch)

    calc_BS!(plank, nuts.T, p, arch)

    update_biomass!(plank, p, ΔT, arch)

    calc_exudation!(plank, p, arch)

    update_CH!(plank, arch)

    ##### probabilities of grazing, mortality, and cell division
    ##### check the probabilities every 10 time steps or 1 hour whichever is shorter
    if t%(ΔT*(min(10,3600÷ΔT))) == 0 
        calc_graz_quadratic!(plank, nuts, p.grz_P, arch)
        calc_MM_mort!(plank, p, arch)
        ##### Bernouli-like distribution
        calc_MM_dvid!(plank, p, arch)
        # plank.dvid .= p.dvid_P .* (1.0 .- isless.(plank.DNA ./ (p.C_DNA * p.Nsuper), 2.0))
        get_probability!(plank, rnd, ΔT, arch)
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end

end