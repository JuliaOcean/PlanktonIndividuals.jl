##### temperature function for photosynthesis
@inline function tempFunc_CF(T, p)
    x = T - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 4.0f-2)
    k = exp(-p.Ea/(8.3145f0*(T+273.15f0))) * regT
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### temperature function for nutrient uptakes
@inline function tempFunc(T, p)
    k = exp(-p.Ea/(8.3145f0*(T+273.15f0)))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, T, Chl, PRO, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * tempFunc_PS(T, p)
    PS  = PCm * (1.0f0 - exp(-αI / max(1.0f-30, PCm) * Chl / max(1.0f-30, PRO))) * PRO
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NP_uptake(NH4, NO3, PO4, T, PRO_Tn, PRO_Tp, PRO_Tfe, pop, p, ac, ΔT)
    VNH4 = p.KcatNH4 * PRO_Tn * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(T, p) * ac
    VNO3 = p.KcatNO3 * PRO_Tn * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(T, p) * ac
    VPO4 = p.KcatPO4 * PRO_Tp * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(T, p) * ac
    VFe  = p.KSAFe * SA * DFe * regQFe * p.Nsuper * tempFunc(T, p) * ac

    return min(VNH4, NH4/ΔT/max(1.0f0,pop)), 
           min(VNO3, NO3/ΔT/max(1.0f0,pop)), 
           min(VPO4, PO4/ΔT/max(1.0f0,pop)),
           min(VFe, DFe/ΔT/max(1.0f0,pop))
end

@kernel function calc_inorganic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(trs.par[i], trs.T[i], plank.Chl[i], plank.PRO[i], p) * plank.ac[i]

    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i], plank.VFe[i] = 
                            calc_NP_uptake(trs.NH4[i], trs.NO3[i], trs.PO4[i], trs.T[i],
                                        plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i],
                                        trs.pop[i], p, plank.ac[i], ΔT)
end
function calc_inorganic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### update C, N, P reserves for the first time of each time step
@kernel function update_quotas_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.CH[i]  += plank.PS[i]   * ΔT
    @inbounds plank.PST[i] += plank.VPO4[i] * ΔT
    @inbounds plank.NST[i] +=(plank.VNH4[i] + plank.VNO3[i]) * ΔT
end
function update_quotas_1!(plank, ΔT, arch)
    kernel! = update_quotas_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 5% of total C acquisition.
@inline function calc_DOC_uptake(DOC, T, CH, PRO, DNA, RNA, Chl, pop, p, ΔT)
    C_tot = total_C_biomass(PRO, DNA, RNA, CH, Chl)
    R_CH = CH / max(1.0f-30, C_tot)
    regQ = shape_func_dec(R_CH, p.CHmax, 1.0f-4)
    VN = p.VDOCmax * regQ * DOC/max(1.0f-30, DOC+p.KsatDOC) * tempFunc(T, p) * PRO
    return min(VN, DOC/ΔT/max(1.0f0,pop))
end
@kernel function calc_organic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(trs.DOC[i], trs.T[i],
                                              plank.CH[i], plank.PRO[i], plank.DNA[i],
                                              plank.RNA[i], plank.Chl[i], 
                                              trs.pop[i], p, ΔT) * plank.ac[i]

    @inbounds plank.VDOC[i] = plank.VDOC[i] * isless(0.05f0, plank.PS[i]/(plank.VDOC[i]+plank.PS[i]))
end
function calc_organic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = plank.PS[i] / max(1.0f-30, plank.PRO[i]) /
                              max(1.0f-30, par[i] * p.α * p.Φ * plank.Chl[i]/max(1.0f-30, plank.PRO[i])) *
                              isless(0.1f0, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, par, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.resp[i] = p.respir * plank.PRO[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end


##### update C, N, P reserves for the second time of each time step
##### respiration first use Carbohydrate, if it's not enough, use protein then.
@kernel function update_quotas_2_kernel!(plank, ΔT, p)
    i = @index(Global)
    @inbounds plank.CH[i]  = plank.CH[i]  + (plank.VDOC[i] - plank.resp[i]) * ΔT
    @inbounds plank.PRO[i] = plank.PRO[i] - max(0.0f0, (0.0f0 - plank.CH[i]))
    @inbounds plank.NST[i] = plank.NST[i] + max(0.0f0, (0.0f0 - plank.CH[i])) * p.R_NC_PRO
    @inbounds plank.CH[i]  = plank.CH[i]  + max(0.0f0, (0.0f0 - plank.CH[i]))
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, T, p)
    i  = @index(Global)
    @inbounds limit_PRO = min(plank.CH[i]/(plank.CH[i] + p.k_sat_pro * p.Nsuper),
                              plank.NST[i]/(plank.NST[i] + p.k_sat_pro * p.Nsuper * p.R_NC_PRO))
    @inbounds limit_DNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_dna * p.Nsuper),
                              plank.NST[i]/(plank.NST[i] + p.k_sat_dna * p.Nsuper * p.R_NC_DNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_dna * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_rna * p.Nsuper),
                              plank.NST[i]/(plank.NST[i] + p.k_sat_rna * p.Nsuper * p.R_NC_RNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_rna * p.Nsuper * p.R_PC_RNA))
    
    @inbounds N_tot = total_N_biomass(plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.DNA[i], plank.RNA[i], plank.NST[i], plank.Chl[i], p)
    @inbounds P_tot = total_P_biomass(plank.DNA[i], plank.RNA[i], plank.PST[i], p)
    @inbounds R_NST = plank.NST[i] / max(1.0f-30, N_tot)
    @inbounds R_PST = plank.PST[i] / max(1.0f-30, P_tot)
    @inbounds regQN = shape_func_dec(R_NST, p.NSTmax, 1.0f-4)
    @inbounds reg_NO3 = shape_fun_dec()
    @inbounds regQP = shape_func_dec(R_PST, p.PSTmax, 1.0f-4)
    @inbounds regQFe = shape_func_dec(QFe, p.qFe_max, 1.0f-4)
    @inbounds Ksat_Fe = QFe / max(1.0f-30, (QFe + p.Ksat_Fe))
     
    @inbounds plank.S_Pr[i] =  plank.PRO_R[i] * p.r_max / p.n_r * limit_PRO * tempFunc(T[i], p) * limit_DNA
    @inbounds plank.S_Pmc[i] =  plank.PRO_Mc[i] * p.r_max / p.n_mC * limit_PRO * tempFunc(T[i], p)
    @inbounds plank.S_Pmn[i] =  plank.PRO_Mn[i] * p.r_max / p.n_mN * limit_PRO * tempFunc(T[i], p) * Ksat_Fe
    @inbounds plank.S_Ptn[i] = plank.PRO_Tn[i] * p.r_max / p.n_tN * limit_PRO * tempFunc(T[i], p) * regQN
    @inbounds plank.S_Ptp[i] =  plank.PRO_Tp[i] * p.r_max / p.n_tPO4 * limit_PRO * tempFunc(T[i], p) * regQP
    @inbounds plank.S_Ptfe[i]= plank.PRO_Tfe[i] * p.r_max / p.n_tFe * limit_PRO * tempFunc(T[i], p) * regQFe 
    @inbounds plank.S_DNA[i] = p.k_dna  * limit_DNA * tempFunc(T[i], p) * isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0f0)
    @inbounds plank.S_RNA[i] = plank.S_Pr[i] * limit_RNA * p.R_C_RNAPr 

function calc_BS!(plank, T, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.PRO_R[i] += ΔT * plank.S_Pr[i]
    @inbounds plank.PRO_Mc[i] += ΔT * plank.S_Pmc[i]
    @inbounds plank.PRO_Mn[i] += ΔT * plank.S_Pmn[i]
    @inbounds plank.PRO_Tn[i] += ΔT * plank.S_Ptn[i]
    @inbounds plank.PRO_Tp[i] += ΔT * plank.S_Ptp[i]
    @inbounds plank.PRO_Tfe[i] += ΔT * plank.S_Ptfe[i]
    @inbounds plank.DNA[i] += ΔT * plank.S_DNA[i]
    @inbounds plank.RNA[i] += ΔT * plank.S_RNA[i]

    @inbounds S_PRO = plank.S_Pr[i] + plank.S_Pmc[i] + plank.S_Pmn[i] + plank.S_Ptn[i] + plank.S_Ptp[i] + plank.S_Ptfe[i]
    @inbounds plank.CH[i]  -= ΔT *(S_PRO + plank.S_DNA[i] + plank.S_RNA[i] + 
                                   plank.S_Pp[i] * plank.ρChl[i])
    @inbounds plank.NST[i] -= ΔT *(S_PRO * p.R_NC_PRO + plank.S_DNA[i] * p.R_NC_DNA + 
                                   plank.S_RNA[i] * p.R_NC_RNA + 
                                   plank.S_Pp[i] * plank.ρChl[i] * 4.0f0 / 55.0f0)
    @inbounds plank.PST[i] -= ΔT *(plank.S_DNA[i] * p.R_PC_DNA + plank.S_RNA[i] * p.R_PC_RNA)
    @inbounds plank.Chl[i] += ΔT * plank.S_Pp[i] * plank.ρChl[i] * 893.49f0 / 55.0f0 # chl unit is mgChl/cell
    @inbounds plank.age[i] += ΔT / 3600.0f0 * plank.ac[i]
end
function update_biomass!(plank, p, ΔT, arch)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### calculate exudation of carbon. Nitrogen and phosphorus will not be exuded for now
@kernel function calc_exudation_kernel!(plank, p)
    i = @index(Global)
    @inbounds tot_C = total_C_biomass(plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i])
    @inbounds plank.exu[i] = max(0.0f0, plank.CH[i] - p.CHmax * tot_C)
end
function calc_exudation!(plank, p, arch)
    kernel! = calc_exudation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

@kernel function update_CH_kernel!(plank)
    i = @index(Global)
    @inbounds plank.CH[i] -= plank.exu[i]
end
function update_CH!(plank, arch)
    kernel! = update_CH_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank)
    return nothing
end
