##### set up the operating array(cuarray) for plankton physiology processes
function phyt_op_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 2)
    op_array = zeros(27, total_num) |> array_type(arch)

    op_array[1:13, :] .= phytos[1:13, :]
    return op_array
end

function grazing!(op_array, consume, arch::Architecture, g::Grids, p)
    grz_array = op_array[:, findall(x -> x == 1.0, op_array[25,:])]
    op_array  = op_array[:, findall(x -> x == 0.0, op_array[25,:])]
    calc_loss!(grz_array, arch, consume.DOC, consume.POC, consume.DON, consume.PON, consume.DOP, consume.POP,
               g, p["grazFracC"], p["grazFracN"], p["grazFracP"], p["_R_NC"], p["R_PC"])
end

function mortality!(op_array, consume, arch::Architecture, g::Grids, p)
    mo_array = op_array[:, findall(x -> x == 1.0, op_array[26,:])]
    op_array = op_array[:, findall(x -> x == 0.0, op_array[26,:])]
    calc_loss!(mo_array, arch, consume.DOC, consume.POC, consume.DON, consume.PON, consume.DOP, consume.POP,
               g, p["mortFracC"], p["mortFracN"], p["mortFracP"], p["_R_NC"], p["R_PC"])
end

function divide!(op_array)
    dvid_array = op_array[:, findall(x -> x == 1.0, op_array[27,:])]
    op_array   = op_array[:, findall(x -> x == 0.0, op_array[27,:])]

    dvid_array[4,:]  .= dvid_array[5,:]  .* 0.45
    dvid_array[5,:]  .= dvid_array[5,:]  .* 0.45
    dvid_array[6,:]  .= dvid_array[6,:]  .* 0.45
    dvid_array[7,:]  .= dvid_array[7,:]  .* 0.5
    dvid_array[8,:]  .= dvid_array[8,:]  .* 0.5
    dvid_array[9,:]  .= dvid_array[9,:]  .* 0.5
    dvid_array[10,:] .= dvid_array[10,:] .* 0.5
    dvid_array[12,:] .= dvid_array[12,:] .+ 1.0
    dvid_array[13,:] .= 1.0

    op_array = hcat(op_array, dvid_array)
    op_array = hcat(op_array, dvid_array)
end


##### update physiological attributes of each individual
function plankton_update!(phytos, consume, diags, arch::Architecture, temp, surf_par, DOC, NH4, NO3, PO4, g, p, ΔT, t)
    chl = zeros(g.Nx, g.Ny, g.Nz) |> array_type(arch)
    par = zeros(g.Nx, g.Ny, g.Nz) |> array_type(arch)
    calc_chla_field!(chl, arch, phytos, g)
    calc_par!(par, arch ,chl, g, surf_par, p["kc"], p["kw"])

    ##### diagnostics for nutrients and par
    diag_t = t ÷ p["diag_freq"]+1
    diags.tr[:,:,:,diag_t,1] += par
    diags.tr[:,:,:,diag_t,2] += NO3
    diags.tr[:,:,:,diag_t,3] += NH4
    diags.tr[:,:,:,diag_t,4] += PO4
    diags.tr[:,:,:,diag_t,5] += DOC

    ##### set up the operating array
    op = phyt_op_array_setup(phytos, arch)

    calc_αI!(op, arch, par, g, p["α"], p["Φ"])
    calc_tempfunc!(op, arch, temp, g, p["TempAe"], p["Tempref"], p["TempCoeff"])

    ##### Carbon uptake
    calc_PS!(op, arch, p["PCmax"], p["PC_b"])
    calc_VDOC!(op, arch, DOC, g, ΔT, p["Cqmax"], p["Cqmin"], p["VDOCmax"], p["VDOC_b"], p["KsatDOC"])

    ##### Nitrogen uptake
    calc_VN!(op, arch, NH4, NO3, g, ΔT,
             p["Nqmax"], p["Nqmin"], p["VNH4max"], p["VNO3max"], p["VN_b"], p["KsatNH4"], p["KsatNO3"], p["R_NC"])

    ##### Phosphorus uptake
    calc_VP!(op, arch, PO4, g, ΔT, p["Pqmax"], p["Pqmin"], p["VPO4max"], p["VP_b"], p["KsatPO4"], p["R_PC"])

    ##### Chla
    calc_ρchl!(op, arch, p["Chl2N"])

    ##### respiration
    calc_respir!(op, arch, p["respir_a"], p["respir_b"])

    ##### update C, N, P quotas
    update_quotas!(op, arch, p["R_NC"], p["R_PC"], ΔT)

    ##### Biosynthesis
    calc_BS!(op, arch, p["k_mtb"], p["b_k_mtb"], p["R_NC"], p["R_PC"])
    updata_biomass!(op, arch, p["R_NC"], p["R_PC"], p["P_Cquota"], p["P_Nsuper"], ΔT)

    calc_consume!(op, arch, consume.DIC, consume.DOC, consume.NH4, consume.NO3, consume.PO4, g, ΔT)

    ##### probabilities of grazing, mortality, and cell division
    if t%300 == 1 # check every 5 mins
        calc_graz!(op, arch, p["Grz_P"], p["Grz_stp"])
        calc_mort!(op, arch, p["mort_reg"], p["mort_P"])
        calc_dvid!(op, arch, p["dvid_type"], p["dvid_stp"], p["dvid_P"], p["divd_reg"], p["P_Cquota"], p["P_Nsuper"])
    end

    ###### diagnostics for each species
    sum_diags!(diags.spcs, op, arch, g, p["diags_inds"], diag_t)

    ##### deal with grazed individual
    grazing!(op, consume, arch, g, p)

    ###### diagnostics of mortality after grazing for each species
    sum_diags_mort!(diags.spcs, op, arch, g, diag_t)

    ##### deal with dead individual
    mortality!(op, consume, arch, g, p)

    ###### diagnostics of cell division after grazing and mortality for each species
    sum_diags_dvid!(diags.spcs, op, arch, g, diag_t)

    ##### deal with divided individual
    divide!(op)

    ##### deal with diagnostics

    phytos = op[1:13, :]
end


