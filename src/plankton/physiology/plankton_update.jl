##### update physiological attributes of each individual
function plankton_update!(plank, rnd, plk, par, arch::Architecture, temp, pop, DOC, NH4, NO3, PO4,
                          g::Grids, p, ΔT, t, plank_num::Int64)
    NO3 = interior(NO3, g)
    NH4 = interior(NH4, g)
    PO4 = interior(PO4, g)
    DOC = interior(DOC, g)

    ##### calculate the total active plankton numbers
    ##### find nutrient, temperature, and par values for each individual
    find_NPT!(plank, Int.(plank[:,13:15]), arch, NH4, NO3, PO4, DOC, par, temp, pop, g,
              p.α, p.Φ, p.TempAe, p.Tempref, p.TempCoeff)

    ##### Carbon uptake
    calc_PS!(plank, arch, p.PCmax, p.PC_b)
    calc_VDOC!(plank, arch, g, ΔT, p.Cqmax, p.Cqmin, p.VDOCmax, p.VDOC_b, p.KsatDOC)

    ##### Nitrogen uptake
    calc_VN!(plank, arch, g, ΔT,
             p.Nqmax, p.Nqmin, p.VNH4max, p.VNO3max, p.VN_b, p.KsatNH4, p.KsatNO3, p.R_NC)

    ##### Phosphorus uptake
    calc_VP!(plank, arch, g, ΔT, p.Pqmax, p.Pqmin, p.VPO4max, p.VP_b, p.KsatPO4, p.R_PC)

    ##### Chla
    calc_ρchl!(plank, arch, p.Chl2N)

    ##### respiration
    calc_respir!(plank, arch, p.respir_a, p.respir_b)

    ##### update C, N, P quotas
    update_quotas!(plank, arch, p.R_NC, p.R_PC, ΔT)

    ##### Biosynthesis
    calc_BS!(plank, arch, p.k_mtb, p.k_mtb_b, p.R_NC, p.R_PC)
    update_biomass!(plank, arch, p.R_NC, p.R_PC, p.Cquota, p.Nsuper, ΔT)

    calc_consume!(plank, Int.(plank[:,13:15]), arch, plk.DIC.data, plk.DOC.data,
                  plk.NH4.data, plk.NO3.data, plk.PO4.data, g, ΔT)

    # ##### probabilities of grazing, mortality, and cell division
    if t%600 == 1
        ##### grazing
        if p.grz_P == 0
            @inbounds plank[1:plank_num,31] .= 0.0
        else
            if p.grz_stp == 0
                calc_graz_quadratic!(plank, arch, p.grz_P)
            else
                calc_graz_linear!(plank, arch, p.grz_P, p.grz_stp)
            end
        end

        ##### mortality
        calc_mort!(plank, arch, p.mort_reg,  p.mort_P)

        ##### cell division
        if p.dvid_type == 1
            calc_dvid_size!(plank, arch, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 2
            calc_dvid_add!(plank, arch, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 3
            calc_dvid_age!(plank, arch, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper)
        elseif p.dvid_type == 4
            calc_dvid_time!(plank, arch, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, t)
        elseif p.dvid_type == 4
            calc_dvid_ts!(plank, arch, p.dvid_type, p.dvid_stp, p.dvid_stp2, p.dvid_P,
                          p.dvid_reg, p.dvid_reg2, p.Cquota, p.Nsuper, t)
        else
            throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
        end

        get_rands!(plank, rnd, arch)
    else
        @inbounds plank[1:plank_num,31] .= 0.0
        @inbounds plank[1:plank_num,32] .= 0.0
        @inbounds plank[1:plank_num,33] .= 0.0
    end
end

