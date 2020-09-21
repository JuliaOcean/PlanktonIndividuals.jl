##### update physiological attributes of each individual
function plankton_update!(plank, plk, par, arch::Architecture, temp, DOC, NH4, NO3, PO4,
                          g::Grids, p, ΔT, t, plank_num::Int64)
    NO3 = interior(NO3, g)
    NH4 = interior(NH4, g)
    PO4 = interior(PO4, g)
    DOC = interior(DOC, g)

    ##### calculate the total active plankton numbers
    ##### find nutrient, temperature, and par values for each individual
    find_NPT!(plank, Int.(plank[:,13:15]), arch, NH4, NO3, PO4, DOC, par, temp, g,
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
                @inbounds plank[1:plank_num,31] .= 1.0/p.grz_P
            else
                @inbounds plank[1:plank_num,31] .=
                    1.0 ./ p.grz_P .* max.(0.15, 1 .- abs.(plank[1:plank_num,3]) ./ p.grz_stp)
            end
        end

        ##### mortality
        @inbounds plank[1:plank_num,32] .= p.mort_P .* (tanh.(6.0 .* (p.mort_reg .- plank[1:plank_num,5])) .+ 1)

        ##### cell division
        calc_dvid!(plank, arch, p.dvid_type, p.dvid_stp, p.dvid_stp2, p.dvid_P,
                   p.dvid_reg, p.dvid_reg2, p.Cquota, p.Nsuper, t)
        get_rands!(plank, arch)
    else
        @inbounds plank[1:plank_num,31] .= 0.0
        @inbounds plank[1:plank_num,32] .= 0.0
        @inbounds plank[1:plank_num,33] .= 0.0
    end
end

