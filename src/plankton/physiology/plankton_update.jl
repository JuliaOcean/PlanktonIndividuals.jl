function grazing!(op_array, plk, arch::Architecture, g::Grids, p)
    grz_array = op_array[findall(x -> x == 1.0, op_array[:,32]), :]
    op_array  = op_array[findall(x -> x == 0.0, op_array[:,32]), :]
    calc_loss!(grz_array, Int.(grz_array[:,14:16]), arch,
               plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p["grazFracC"], p["grazFracN"], p["grazFracP"], p["_R_NC"], p["R_PC"])
end

function mortality!(op_array, plk, arch::Architecture, g::Grids, p)
    mo_array = op_array[findall(x -> x == 1.0, op_array[:,33]), :]
    op_array = op_array[findall(x -> x == 0.0, op_array[:,33]), :]
    calc_loss!(mo_array, Int.(mo_array[:,14:16]), arch,
               plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p["mortFracC"], p["mortFracN"], p["mortFracP"], p["_R_NC"], p["R_PC"])
end

function divide!(op_array)
    dvid_array = op_array[findall(x -> x == 1.0, op_array[:,31]), :]
    op_array   = op_array[findall(x -> x == 0.0, op_array[:,31]), :]

    dvid_array[:,4]  .= dvid_array[:,5]  .* 0.45
    dvid_array[:,5]  .= dvid_array[:,5]  .* 0.45
    dvid_array[:,6]  .= dvid_array[:,6]  .* 0.45
    dvid_array[:,7]  .= dvid_array[:,7]  .* 0.5
    dvid_array[:,8]  .= dvid_array[:,8]  .* 0.5
    dvid_array[:,9]  .= dvid_array[:,9]  .* 0.5
    dvid_array[:,10] .= dvid_array[:,10] .* 0.5
    dvid_array[:,12] .= dvid_array[:,12] .+ 1.0
    dvid_array[:,13] .= 1.0

    op_array = hcat(op_array, dvid_array)
    op_array = hcat(op_array, dvid_array)
end


##### update physiological attributes of each individual
function plankton_update!(phytos, ope, plk, par, chl, diags, arch::Architecture,
                          temp, surf_par, DOC, NH4, NO3, PO4, g::Grids, p, ΔT, t)
    NO3 = interior(NO3.data, g)
    NH4 = interior(NH4.data, g)
    PO4 = interior(PO4.data, g)
    DOC = interior(DOC.data, g)

    ##### copy phytos to the operating array
    update_phy_ope!(phytos, ope, arch)

    ##### find inds
    find_inds!(ope, arch, g, 13, 0)

    ##### calculate chl and par fields
    acc_chla_field!(chl, ope, Int.(ope[:,14:16]), arch)
    calc_par!(par, arch, chl, g, surf_par, p["kc"], p["kw"])

    # ##### diagnostics for nutrients and par
    # diag_t = t ÷ p["diag_freq"] + 1
    # diags.tr[:,:,:,diag_t,1] += par
    # diags.tr[:,:,:,diag_t,2] += NO3
    # diags.tr[:,:,:,diag_t,3] += NH4
    # diags.tr[:,:,:,diag_t,4] += PO4
    # diags.tr[:,:,:,diag_t,5] += DOC

    ##### find nutrient, temperature, and par values for each individual
    find_NPT!(ope, Int.(ope[:,14:16]), Int.(ope[:,11]), arch, NH4, NO3, PO4, DOC, par, temp, g,
              p["α"], p["Φ"], p["TempAe"], p["Tempref"], p["TempCoeff"])

    ##### Carbon uptake
    calc_PS!(ope, Int.(ope[:,11]), arch, p["PCmax"], p["PC_b"])
    calc_VDOC!(ope, Int.(ope[:,11]), arch, g, ΔT, p["Cqmax"], p["Cqmin"], p["VDOCmax"], p["VDOC_b"], p["KsatDOC"])

    ##### Nitrogen uptake
    calc_VN!(ope, Int.(ope[:,11]), arch, g, ΔT,
             p["Nqmax"], p["Nqmin"], p["VNH4max"], p["VNO3max"], p["VN_b"], p["KsatNH4"], p["KsatNO3"], p["R_NC"])

    ##### Phosphorus uptake
    calc_VP!(ope, Int.(ope[:,11]), arch, g, ΔT, p["Pqmax"], p["Pqmin"], p["VPmax"], p["VP_b"], p["KsatP"], p["R_PC"])

    ##### Chla
    calc_ρchl!(ope, arch, p["Chl2N"])

    ##### respiration
    calc_respir!(ope, Int.(ope[:,11]), arch, p["respir_a"], p["respir_b"])

    ##### update C, N, P quotas
    update_quotas!(ope, arch, p["R_NC"], p["R_PC"], ΔT)

    ##### Biosynthesis
    calc_BS!(ope, Int.(ope[:,11]), arch, p["k_mtb"], p["b_k_mtb"], p["R_NC"], p["R_PC"])
    update_biomass!(ope, Int.(ope[:,11]), arch, p["R_NC"], p["R_PC"], p["P_Cquota"], p["P_Nsuper"], ΔT)

    calc_consume!(ope, Int.(ope[:,14:16]), arch, plk.DIC.data, plk.DOC.data,
                  plk.NH4.data, plk.NO3.data, plk.PO4.data, g, ΔT)

    # ##### probabilities of grazing, mortality, and cell division
    # probabilities!(ope, Int.(ope[:,11]), arch, p["Grz_P"], p["Grz_stp"], p["mort_reg"], p["mort_P"],
    #                p["dvid_type"], p["dvid_stp"], p["dvid_stp2"], p["dvid_P"], p["dvid_reg"],
    #                p["dvid_reg2"], p["P_Cquota"], p["P_Nsuper"], t)

    # ###### diagnostics for each species
    # sum_diags!(diags.spcs, ope, Int.(ope[:,14:16]), Int.(ope[:,11]), arch, g, p["diag_inds"], diag_t)

    # ##### deal with grazed individual
    # grazing!(ope, plk, arch, g, p)

    # ###### diagnostics of mortality after grazing for each species
    # sum_diags_mort!(diags.spcs, ope, Int.(ope[:,14:16]), Int.(ope[:,11]), arch, g, diag_t)

    # ##### deal with dead individual
    # mortality!(ope, plk, arch, g, p)

    # ###### diagnostics of cell division after grazing and mortality for each species
    # sum_diags_dvid!(diags.spcs, ope, Int.(ope[:,14:16]), Int.(ope[:,11]), arch, g, diag_t)

    # ##### deal with divided individual
    # divide!(ope)

    phytos = ope[:, 1:13]
end


