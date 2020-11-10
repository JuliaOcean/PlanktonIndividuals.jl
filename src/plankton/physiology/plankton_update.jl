##### update physiological attributes of each individual
function plankton_update!(plank, nuts, coord, tmp, rnd, cts, par, arch::Architecture,
                          temp, pop, nut, g::Grids, p, ΔT, t, num::Int64)
    ##### calculate the total active plankton numbers
    ##### find nutrient, temperature, and par values for each individual
    find_NPT!(nuts, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g,
              nut.NH4.data, nut.NO3.data, nut.PO4.data, nut.DOC.data,
              par, temp, pop, p.α, p.Φ, p.TempAe, p.Tempref, p.TempCoeff)

    # ##### Carbon uptake
    # calc_PS!(plank, tmp, p.PCmax, p.PC_b, num)
    # calc_VDOC!(plank, tmp, g, ΔT, p.Cqmax, p.Cqmin, p.VDOCmax, p.VDOC_b, p.KsatDOC, num)

    # ##### Nitrogen uptake
    # calc_VN!(plank, tmp, g, ΔT, p.Nqmax, p.Nqmin,
    #          p.VNH4max, p.VNO3max, p.VN_b, p.KsatNH4, p.KsatNO3, p.R_NC, num)

    # ##### Phosphorus uptake
    # calc_VP!(plank, tmp, g, ΔT, p.Pqmax, p.Pqmin, p.VPO4max, p.VP_b, p.KsatPO4, p.R_PC, num)

    # ##### Chla
    # calc_ρchl!(plank, tmp, p.Chl2N, num)

    # ##### respiration
    # calc_respir!(plank, p.respir_a, p.respir_b, num)

    # ##### update C, N, P quotas
    # update_quotas!(plank, tmp, p.R_NC, p.R_PC, ΔT, num)

    # ##### Biosynthesis
    # calc_BS!(plank, tmp, p.k_mtb, p.k_mtb_b, p.R_NC, p.R_PC, num)
    # update_biomass!(plank, p.R_NC, p.R_PC, p.Cquota, p.Nsuper, ΔT, num)

    # calc_consume!(cts, plank, Int.(plank[:,13:15]), arch, g, ΔT)

    # # ##### probabilities of grazing, mortality, and cell division
    # if t%600 == 1
    #     ##### grazing
    #     if p.grz_P == 0
    #         @inbounds plank[1:num,31] .= 0.0
    #     else
    #         if p.grz_stp == 0
    #             calc_graz_quadratic!(plank, p.grz_P, num)
    #         else
    #             calc_graz_linear!(plank, p.grz_P, p.grz_stp, num)
    #         end
    #     end

    #     ##### mortality
    #     calc_mort!(plank, p.mort_reg,  p.mort_P, num)

    #     ##### cell division
    #     if p.dvid_type == 1
    #         calc_dvid_size!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, num)
    #     elseif p.dvid_type == 2
    #         calc_dvid_add!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, num)
    #     elseif p.dvid_type == 3
    #         calc_dvid_age!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, num)
    #     elseif p.dvid_type == 4
    #         calc_dvid_time!(plank, p.dvid_stp, p.dvid_P, p.dvid_reg, p.Cquota, p.Nsuper, t, num)
    #     elseif p.dvid_type == 4
    #         calc_dvid_ts!(plank, p.dvid_type, p.dvid_stp, p.dvid_stp2, p.dvid_P,
    #                       p.dvid_reg, p.dvid_reg2, p.Cquota, p.Nsuper, t, num)
    #     else
    #         throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
    #     end

    #     get_rands!(plank, rnd, num)
    # else
    #     @inbounds plank[1:num,31] .= 0.0
    #     @inbounds plank[1:num,32] .= 0.0
    #     @inbounds plank[1:num,33] .= 0.0
    # end
end

