function grazing!(plank, plk, arch::Architecture, g::Grids, p)
    grz_array = plank[findall(x -> x == 1.0, plank[:,31]), :]
    calc_loss!(grz_array, Int.(grz_array[:,13:15]), arch,
               plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.grazFracC, p.grazFracN, p.grazFracP, p.R_NC, p.R_PC)
    plank = plank[findall(x -> x == 0.0, plank[:,31]), :]
end

function mortality!(plank, plk, arch::Architecture, g::Grids, p)
    mo_array = plank[findall(x -> x == 1.0, plank[:,32]), :]

    calc_loss!(mo_array, Int.(mo_array[:,13:15]), arch,
               plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.mortFracC, p.mortFracN, p.mortFracP, p.R_NC, p.R_PC)

    plank = plank[findall(x -> x == 0.0, plank[:,32]), :]
end

function divide!(plank)
    dvid_array = plank[findall(x -> x == 1.0, plank[:,33]), :]
    plank   = plank[findall(x -> x == 0.0, plank[:,33]), :]

    dvid_array[:,4]  .= dvid_array[:,5]  .* 0.45
    dvid_array[:,5]  .= dvid_array[:,5]  .* 0.45
    dvid_array[:,6]  .= dvid_array[:,6]  .* 0.45
    dvid_array[:,7]  .= dvid_array[:,7]  .* 0.5
    dvid_array[:,8]  .= dvid_array[:,8]  .* 0.5
    dvid_array[:,9]  .= dvid_array[:,9]  .* 0.5
    dvid_array[:,10] .= dvid_array[:,10] .* 0.5
    dvid_array[:,11] .= dvid_array[:,11] .+ 1.0
    dvid_array[:,12] .= 1.0

    plank = vcat(plank, dvid_array)
    plank = vcat(plank, dvid_array)
end


##### update physiological attributes of each individual
function plankton_update!(plank, plk, par, arch::Architecture,
                          temp, DOC, NH4, NO3, PO4, g::Grids, p, ΔT, t)
    NO3 = interior(NO3, g)
    NH4 = interior(NH4, g)
    PO4 = interior(PO4, g)
    DOC = interior(DOC, g)

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
            @inbounds plank[:,31] .= 0.0
        else
            if p.grz_stp == 0
                @inbounds plank[:,31] .= 1.0/p.grz_P
            else
                @inbounds plank[:,31] .= 1.0 ./ p.grz_P .* max.(0.15, 1 .- abs.(plank[:,3]) ./ p.grz_stp)
            end
        end

        ##### mortality
        @inbounds plank[:,32] .= p.mort_P .* (tanh.(6.0 .* (p.mort_reg .- plank[:,5])) .+ 1)

        ##### cell division
        calc_dvid!(plank, arch, p.dvid_type, p.dvid_stp, p.dvid_stp2, p.dvid_P,
                   p.dvid_reg, p.dvid_reg2, p.Cquota, p.Nsuper, t)
        get_rands!(plank, arch)
    else
        @inbounds plank[:,31] .= 0.0
        @inbounds plank[:,32] .= 0.0
        @inbounds plank[:,33] .= 0.0
    end
end


