##### update physiological attributes of each individual
function plankton_update!(plank, nuts, proc, rnd, par, pop,
                          temp, nut, g::Grids, p, ΔT, t)
    ##### calculate the total active plankton numbers
    ##### find nutrient, temperature, and par values for each individual
    find_NPT!(nuts, Int.(plank.xi), Int.(plank.yi), Int.(plank.zi), plank.ac, g,
              nut.NH4.data, nut.NO3.data, nut.PO4.data, nut.DOC.data,
              par, temp, pop, p.α, p.Φ, p.TempAe, p.Tempref, p.TempCoeff)

    ##### Carbon uptake
    calc_PS!(plank, proc, nuts, p.PCmax, p.PC_b)

    ##### Nitrogen uptake
    calc_VN!(plank, proc, nuts, g, ΔT, p.Nqmax, p.Nqmin,
             p.VNH4max, p.VNO3max, p.VN_b, p.KsatNH4, p.KsatNO3, p.R_NC)

    ##### Phosphorus uptake
    calc_VP!(plank, proc, nuts, g, ΔT, p.Pqmax, p.Pqmin, p.VPO4max, p.VP_b, p.KsatPO4, p.R_PC)

    ##### Chla
    calc_ρchl!(plank, proc, nuts, p.Chl2N)

    ##### respiration
    calc_respir!(plank, proc, nuts, p.respir_a, p.respir_b)

    ##### update C, N, P quotas
    update1_quotas!(plank, proc, ΔT)

    ##### Organic carbon uptake if necessary
    calc_VDOC!(plank, proc, nuts, g, ΔT, p.Cqmax, p.Cqmin, p.VDOCmax, p.VDOC_b, p.KsatDOC)

    ##### update C, N, P quotas
    update2_quotas!(plank, proc, p.R_NC, p.R_PC, ΔT)

    ##### Biosynthesis
    calc_BS!(plank, proc, p.k_mtb, p.k_mtb_b, p.R_NC, p.R_PC)
    update_biomass!(plank, proc, p.R_NC, p.R_PC, p.Cquota, p.Nsuper, ΔT)

    ##### probabilities of grazing, mortality, and cell division
    if t%600 == 1
        ##### grazing
        if p.grz_P == 0
            @inbounds plank.graz .= 0.0
        else
            if p.grz_stp == 0
                calc_graz_quadratic!(nuts, plank, p.grz_P)
            else
                calc_graz_linear!(plank, p.grz_P, p.grz_stp)
            end
        end

        ##### mortality
        calc_mort!(plank, p.mort_reg,  p.mort_P)

        ##### cell division
        if p.dvid_type == 1
            calc_dvid_size!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 2
            calc_dvid_add!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 3
            calc_dvid_age!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 4
            calc_dvid_time!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, t)
        elseif p.dvid_type == 5
            calc_dvid_ts!(plank, p.dvid_type, p.dvid_stp, p.dvid_stp2, p.dvid_P,
                          p.dvid_reg, p.dvid_reg2, p.Cquota, p.Nsuper, t)
        else
            throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
        end

        get_rands!(plank, rnd)
    else
        @inbounds plank.graz .= 0.0
        @inbounds plank.mort .= 0.0
        @inbounds plank.dvid .= 0.0
    end
end

