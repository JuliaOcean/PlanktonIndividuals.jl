##### temperature function for carbon fixation
@inline function tempFunc_CF(T, p)
    x = T - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 4.0f-2)
    k = exp(-p.Ea / (8.3145f0 * (T + 273.15f0))) * regT
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea / (8.3145f0 * (p.Topt + 273.15f0)))
    return min(1.0f0, k / OGT_rate)
end

##### temperature function for nutrient uptakes
@inline function tempFunc(T, p)
    k = exp(-p.Ea / (8.3145f0 * (T + 273.15f0)))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea / (8.3145f0 * (p.Topt + 273.15f0)))
    return min(1.0f0, k / OGT_rate)
end

##### calculate photosynthesis rate (mmolATP/individual/second)
@kernel function calc_PS_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds αI = trs.par[i] * p.α * plank.Chl[i] / max(1.0f-30, plank.PRO_PS[i])
    @inbounds plank.PS[i] = p.PCmax * (1.0f0 - exp(-αI)) * plank.PRO_PS[i] * plank.ac[i]
    @inbounds plank.exE_PS[i] = copy(plank.PS[i])
end
function calc_PS!(plank, trs, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p, ΔT)
    i = @index(Global)
    @inbounds plank.RS[i] = p.KcatRS * plank.PRO_RS[i] * tempFunc(T[i], p) * plank.ac[i]
    @inbounds plank.RS[i] = min(plank.RS[i], plank.CH[i] / ΔT)
    @inbounds plank.ERS[i]= plank.RS[i] * p.e_RS * plank.ac[i]
    @inbounds plank.exE_RS[i] = copy(plank.ERS[i])
end
function calc_respiration!(plank, T, p, ΔT, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p, ΔT)
    return nothing
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@kernel function calc_inorganic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])

    @inbounds R_qNH4 = plank.qNH4[i] / max(1.0f-30, C_tot)
    @inbounds R_qNO3 = plank.qNO3[i] / max(1.0f-30, C_tot)
    @inbounds R_PST  = plank.PST[i]  / max(1.0f-30, C_tot)
    @inbounds R_qFe  = plank.qFe[i]  / max(1.0f-30, C_tot)

    @inbounds regQNH4 = shape_func_dec(R_qNH4, p.qNH4max, 1.0f-4)
    @inbounds regQNO3 = shape_func_dec(R_qNO3, p.qNO3max, 1.0f-4)
    @inbounds regQP   = shape_func_dec(R_PST,  p.PSTmax,  1.0f-4)
    @inbounds regQFe  = shape_func_dec(R_qFe,  p.qFemax,  1.0f-4)

    @inbounds MM_NH4 = trs.NH4[i] / max(1.0f-30, trs.NH4[i] + p.KsatNH4)
    @inbounds MM_NO3 = trs.NO3[i] / max(1.0f-30, trs.NO3[i] + p.KsatNO3)
    @inbounds MM_PO4 = trs.PO4[i] / max(1.0f-30, trs.PO4[i] + p.KsatPO4)
    @inbounds MM_DFe = trs.DFe[i] / max(1.0f-30, trs.DFe[i] + p.KsatFe)

    @inbounds plank.VNH4[i] = p.KcatTNH4 * plank.PRO_TN[i]  * MM_NH4 * regQNH4
    @inbounds plank.VNO3[i] = p.KcatTNO3 * plank.PRO_TN[i]  * MM_NO3 * regQNO3
    @inbounds plank.VPO4[i] = p.KcatTPO4 * plank.PRO_TP[i]  * MM_PO4 * regQP
    @inbounds plank.VFe[i]  = p.KcatTFe  * plank.PRO_TFe[i] * MM_DFe * regQFe

    @inbounds plank.VNH4[i] *= tempFunc(trs.T[i], p) * plank.ac[i]  #mmolN/individual/second
    @inbounds plank.VNO3[i] *= tempFunc(trs.T[i], p) * plank.ac[i]  #mmolN/individual/second
    @inbounds plank.VPO4[i] *= tempFunc(trs.T[i], p) * plank.ac[i]  #mmolP/individual/second
    @inbounds plank.VFe[i]  *= tempFunc(trs.T[i], p) * plank.ac[i]  #mmolFe/individual/second

    @inbounds plank.VNH4[i] = min(plank.VNH4[i], trs.NH4[i]/ΔT/max(1.0f0, trs.pop[i]))
    @inbounds plank.VNO3[i] = min(plank.VNO3[i], trs.NO3[i]/ΔT/max(1.0f0, trs.pop[i]))
    @inbounds plank.VPO4[i] = min(plank.VPO4[i], trs.PO4[i]/ΔT/max(1.0f0, trs.pop[i]))
    @inbounds plank.VFe[i]  = min(plank.VFe[i],  trs.DFe[i]/ΔT/max(1.0f0, trs.pop[i]))
    
    @inbounds plank.EVNO3[i] = plank.VNO3[i] * p.e_TNO3
    @inbounds plank.EVPO4[i] = plank.VPO4[i] * p.e_TPO4
    @inbounds plank.EVFe[i]  = plank.VFe[i]  * p.e_TFe 
end
function calc_inorganic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate energy consumption of nutrient uptakes (mmolATP/individual/second)
@kernel function calc_uptake_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds exE_RS = max(plank.exE_RS[i]- max(p.e_min - plank.exE_PS[i], 0.0f0), 0.0f0)
    @inbounds exE_PS = max(plank.exE_PS[i]- p.e_min, 0.0f0)
    @inbounds Esupply = plank.exE_RS[i]
    @inbounds Edemand = plank.EVNO3[i] + plank.EVPO4[i] + plank.EVFe[i]
    @inbounds Eused = min(Esupply, Edemand)

    @inbounds plank.EVNO3[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.EVPO4[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.EVFe[i]  *= Eused / max(1.0f-30, Edemand)

    @inbounds plank.exE_RS[i] -= Eused

    @inbounds plank.VNO3[i] = plank.EVNO3[i] / p.e_TNO3
    @inbounds plank.VPO4[i] = plank.EVPO4[i] / p.e_TPO4
    @inbounds plank.VFe[i]  = plank.EVFe[i]  / p.e_TFe
end

function calc_uptake_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_uptake_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P reserves for the first time of each time step
@kernel function update_quotas_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.PST[i]  += plank.VPO4[i] * ΔT
    @inbounds plank.qNH4[i] += plank.VNH4[i] * ΔT
    @inbounds plank.qNO3[i] += plank.VNO3[i] * ΔT
    @inbounds plank.qFe[i]  += plank.VFe[i]  * ΔT
    @inbounds plank.CH[i]   -= plank.RS[i]   * ΔT
end
function update_quotas_1!(plank, ΔT, arch)
    kernel! = update_quotas_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate potential carbon fixation rate (mmolC/individual/second)
##### and energy consumption rate of carbon fixation (mmolATP/individual/second)
@kernel function calc_carbon_fixation_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])

    @inbounds regQC = shape_func_dec(plank.CH[i] / max(1.0f-30, C_tot), p.CHmax, 1.0f-4)

    @inbounds plank.CF[i]  = p.KcatCF * plank.PRO_MC[i] * regQC
    @inbounds plank.CF[i] *= tempFunc_CF(T[i], p) * plank.ac[i]

    @inbounds plank.ECF[i] = plank.CF[i] * p.e_CF
end
function calc_carbon_fixation!(plank, T, p, arch)
    kernel! = calc_carbon_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate potential nitrate reduction (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@kernel function calc_NO3_reduction_kernel!(plank, T, p, ΔT)
    i = @index(Global)
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])

    @inbounds Qn = plank.qNH4[i]/max(1.0f-30, C_tot)
    @inbounds reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    @inbounds Ksat = plank.qNO3[i] / max(1.0f-30, (plank.qNO3[i] + p.KsatNR))

    @inbounds plank.NR[i]  = p.KcatNR * plank.PRO_MN[i] * reg * Ksat 
    @inbounds plank.NR[i] *= tempFunc(T[i], p) * plank.ac[i] * p.is_nr
    @inbounds plank.NR[i]  = min(plank.NR[i], plank.qNO3[i]/ΔT) # double check qNO3 are not over consumed

    @inbounds plank.ENR[i] = plank.NR[i] * p.e_NR
end
function calc_NO3_reduction!(plank, T, p, ΔT, arch::Architecture)
    kernel! = calc_NO3_reduction_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p, ΔT)
    return nothing
end

##### calculate potential nitrogen fixation (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@kernel function calc_nitrogen_fixation_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])

    @inbounds Qn = plank.qNH4[i]/max(1.0f-30, C_tot)
    @inbounds reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    
    @inbounds plank.NF[i]  = p.KcatNF * plank.PRO_MN[i] * reg
    @inbounds plank.NF[i] *= tempFunc(T[i], p) * plank.ac[i] * (p.is_croc + p.is_tric)

    @inbounds plank.ENF[i] = plank.NF[i] * p.e_NF                               
end
function calc_nitrogen_fixation!(plank, T, p, arch::Architecture)
    kernel! = calc_nitrogen_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate energy consumption of carbon and nitrigen metabolism (mmolATP/individual/second)
@kernel function calc_CN_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.ECF[i] = min(plank.ECF[i], plank.exE_PS[i])
    @inbounds plank.ENR[i] = min(plank.ENR[i], plank.exE_PS[i] + plank.exE_RS[i] - plank.ECF[i])
    @inbounds plank.ENF[i] = min(plank.ENF[i], plank.exE_PS[i] + plank.exE_RS[i] - plank.ECF[i])

    @inbounds plank.exE_RS[i] -= max(0.0f0, plank.ECF[i] + plank.ENR[i] + plank.ENF[i] - plank.exE_PS[i])
    @inbounds plank.exE_RS[i] = max(0.0f0, plank.exE_RS[i])
    @inbounds plank.exE_PS[i] -= min(plank.exE_PS[i], plank.ECF[i] + plank.ENR[i] + plank.ENF[i])
       
    @inbounds plank.CF[i] = plank.ECF[i] / p.e_CF
    @inbounds plank.NR[i] = plank.ENR[i] / p.e_NR
    @inbounds plank.NF[i] = plank.ENF[i] / p.e_NF 
end
function calc_CN_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_CN_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 5% of total C acquisition.
@kernel function calc_organic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])
    
    @inbounds R_CH = plank.CH[i] / max(1.0f-30, C_tot)
    @inbounds regQ = shape_func_dec(R_CH, p.CHmax, 1.0f-4)
    @inbounds MM_DOC = trs.DOC[i] / max(1.0f-30, trs.DOC[i] + p.KsatDOC)
    @inbounds plank.VDOC[i]  = p.VDOCmax * regQ * MM_DOC * C_tot
    @inbounds plank.VDOC[i] *= tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.VDOC[i] *= isless(0.05f0, plank.CF[i] / max(1.0f-30, plank.VDOC[i] + plank.CF[i]))
    @inbounds plank.VDOC[i]  = min(plank.VDOC[i], trs.DOC[i]/ΔT/max(1.0f0, trs.pop[i]))
end
function calc_organic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = p.Chl2N / max(1.0f-30, par[i] * p.α * plank.Chl[i] / 
                              max(1.0f-30, plank.PRO_PS[i])) *
                              isless(1.0f-1, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, par, p)
    return nothing
end

##### calculate protein and chla degration under N stravation (mmolC/individual/second)
@kernel function calc_degradation_kernel!(plank, p, trs)
    i = @index(Global)
    @inbounds reg_PRB = shape_func_dec(plank.qNH4[i], p.qNH4max, 1.0f-4, pow = 1.0f0)
    @inbounds reg_PPS = shape_func_dec(plank.qNH4[i], p.qNH4max, 1.0f-4, pow = 4.0f0)
    @inbounds reg_PMC = shape_func_dec(plank.qNH4[i], p.qNH4max, 1.0f-4, pow = 2.0f0)
    @inbounds reg_Chl = shape_func_dec(plank.qNH4[i], p.qNH4max, 1.0f-4, pow = 4.0f0)

    @inbounds plank.DP_RB[i] = p.k_degRB * reg_PRB * (plank.PRO_RB[i] - p.PRO_RBmin)
    @inbounds plank.DP_PS[i] = p.k_degPS * reg_PPS * (plank.PRO_PS[i] - p.PRO_PSmin)
    @inbounds plank.DP_MC[i] = p.k_degMC * reg_PMC * (plank.PRO_MC[i] - p.PRO_MCmin)
    @inbounds plank.DChl[i]  = p.k_degChl* reg_Chl * (plank.Chl[i] - p.Chlmin)

    @inbounds plank.DP_RB[i] *= tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.DP_PS[i] *= tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.DP_MC[i] *= tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.DChl[i]  *= tempFunc(trs.T[i], p) * plank.ac[i]
end 
function calc_degradation!(plank, p, trs, arch)
    kernel! = calc_degradation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, trs)
    return nothing
end

##### update C, N, P reserves for the second time of each time step
##### respiration first use Carbohydrate, if it's not enough, use protein then.
@kernel function update_quotas_2_kernel!(plank, ΔT, p)
    i = @index(Global)
    @inbounds plank.CH[i] += ΔT * (plank.CF[i] + plank.VDOC[i])
    @inbounds plank.CH[i] += ΔT * (plank.DP_RB[i] + plank.DP_PS[i] + plank.DP_MC[i])
    @inbounds plank.CH[i] += ΔT * plank.DChl[i] / 893.49f0 * 55.0f0

    @inbounds plank.qNH4[i] += ΔT * (plank.NR[i] + plank.NF[i])
    @inbounds plank.qNH4[i] += ΔT * (plank.DP_RB[i] + plank.DP_PS[i] + plank.DP_MC[i]) * p.R_NC_PRO
    @inbounds plank.qNH4[i] += ΔT * plank.DChl[i] / 893.49f0 * 4.0f0

    @inbounds plank.qNO3[i]   -= ΔT * plank.NR[i]
    @inbounds plank.PRO_RB[i] -= ΔT * plank.DP_RB[i]
    @inbounds plank.PRO_PS[i] -= ΔT * plank.DP_PS[i]
    @inbounds plank.PRO_MC[i] -= ΔT * plank.DP_MC[i]
    @inbounds plank.Chl[i]    -= ΔT * plank.DChl[i]
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, trs, p)
    i  = @index(Global)
    @inbounds NST = plank.qNH4[i] + plank.qNO3[i]
    @inbounds limit_PRO = min(plank.CH[i]/(plank.CH[i] + p.k_sat_PRO * p.Nsuper),
                              NST/(NST + p.k_sat_PRO * p.Nsuper * p.R_NC_PRO))
    @inbounds limit_DNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_DNA * p.Nsuper),
                              NST/(NST + p.k_sat_DNA * p.Nsuper * p.R_NC_DNA),
                              plank.PST[i]/(plank.PST[i] + p.k_sat_DNA * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_RNA * p.Nsuper),
                              NST/(NST + p.k_sat_RNA * p.Nsuper * p.R_NC_RNA),
                              plank.PST[i]/(plank.PST[i] + p.k_sat_RNA * p.Nsuper * p.R_PC_RNA))
    
    @inbounds C_tot = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])
    @inbounds QFe  = plank.qFe[i] / max(1.0f-30, C_tot)
    @inbounds Ksat_Fe = QFe / max(1.0f-30, (QFe + p.KFe_em))
    @inbounds regI = trs.par[i] * p.α * plank.Chl[i] / max(1.0f-30, plank.PRO_PS[i])
    @inbounds regI = 1.0f0 / max(1.0f-30, regI)
    @inbounds lim_TN =  shape_func_dec(trs.NO3[i] + trs.NH4[i], p.TN_max, 1.0f-4)
    @inbounds lim_TP =  shape_func_dec(trs.PO4[i], p.TP_max, 1.0f-4)
    @inbounds lim_TFe = shape_func_dec(trs.DFe[i], p.TFe_max, 1.0f-4)
     
    @inbounds plank.SP_RB[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_RB   * limit_RNA
    @inbounds plank.SP_MC[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_MC
    @inbounds plank.SP_MN[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_MN   * Ksat_Fe
    @inbounds plank.SP_TN[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_TN   * lim_TN
    @inbounds plank.SP_TP[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_TPO4 * lim_TP
    @inbounds plank.SP_TFe[i]= plank.PRO_RB[i]  * p.KcatRB * p.β_TFe  * lim_TFe
    @inbounds plank.SP_RS[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_RS   * Ksat_Fe
    @inbounds plank.SP_PS[i] = plank.PRO_RB[i]  * p.KcatRB * p.β_PS   * regI * Ksat_Fe

    @inbounds plank.SP_RB[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_MC[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_MN[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_TN[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_TP[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_TFe[i]*= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_RS[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.SP_PS[i] *= limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    
    @inbounds plank.SDNA[i] = p.k_DNA  * limit_DNA * plank.ac[i] * tempFunc(trs.T[i], p) 
    @inbounds plank.SDNA[i] *= isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0f0)
    @inbounds plank.SRNA[i] = plank.SP_RB[i] * p.R_C_RNAPRB * plank.ac[i] 

    @inbounds plank.ESP_RB[i] = plank.SP_RB[i] * p.e_SP
    @inbounds plank.ESP_MC[i] = plank.SP_MC[i] * p.e_SP
    @inbounds plank.ESP_MN[i] = plank.SP_MN[i] * p.e_SP
    @inbounds plank.ESP_TN[i] = plank.SP_TN[i] * p.e_SP
    @inbounds plank.ESP_TP[i] = plank.SP_TP[i] * p.e_SP
    @inbounds plank.ESP_TFe[i]= plank.SP_TFe[i]* p.e_SP
    @inbounds plank.ESP_RS[i] = plank.SP_RS[i] * p.e_SP
    @inbounds plank.ESP_PS[i] = plank.SP_PS[i] * p.e_SP
    @inbounds plank.EDNA[i]  = plank.SDNA[i] * p.e_DNA
    @inbounds plank.ERNA[i]  = plank.SRNA[i] * p.e_RNA  
end
function calc_BS!(plank, trs, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate energy consumption of protein, DNA, RNA synthesis (mmolATP/individual/second)
@kernel function calc_BS_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds Esupply = plank.exE_PS[i] + plank.exE_RS[i]
    @inbounds Edemand = plank.ESP_RB[i] + plank.ESP_MC[i] + plank.ESP_MN[i] + plank.ESP_TN[i] +
                        plank.ESP_TP[i] + plank.ESP_TFe[i]+ plank.ESP_RS[i] + plank.ESP_PS[i]
    @inbounds Eused = min(Esupply, Edemand)
    
    @inbounds plank.ESP_RB[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_MC[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_MN[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_TN[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_TP[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_TFe[i]*= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_RS[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ESP_PS[i] *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.EDNA[i]   *= Eused / max(1.0f-30, Edemand)
    @inbounds plank.ERNA[i]   *= Eused / max(1.0f-30, Edemand)

    @inbounds plank.exE_RS[i] -= min(plank.exE_RS[i], max(Eused - plank.exE_PS[i], 0.0f0))
    @inbounds plank.exE_PS[i] -= min(plank.exE_PS[i], Eused)
    
    @inbounds plank.SP_RB[i] = plank.ESP_RB[i]  / p.e_SP
    @inbounds plank.SP_MC[i] = plank.ESP_MC[i]  / p.e_SP
    @inbounds plank.SP_MN[i] = plank.ESP_MN[i]  / p.e_SP
    @inbounds plank.SP_TN[i] = plank.ESP_TN[i]  / p.e_SP
    @inbounds plank.SP_TP[i] = plank.ESP_TP[i]  / p.e_SP
    @inbounds plank.SP_TFe[i]= plank.ESP_TFe[i] / p.e_SP
    @inbounds plank.SP_RS[i] = plank.ESP_RS[i]  / p.e_SP
    @inbounds plank.SP_PS[i] = plank.ESP_PS[i]  / p.e_SP
    @inbounds plank.SDNA[i]  = plank.EDNA[i]    / p.e_DNA
    @inbounds plank.SRNA[i]  = plank.ERNA[i]    / p.e_RNA
end
function calc_BS_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_BS_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.PRO_RB[i]  += ΔT * plank.SP_RB[i]
    @inbounds plank.PRO_MC[i]  += ΔT * plank.SP_MC[i]
    @inbounds plank.PRO_MN[i]  += ΔT * plank.SP_MN[i]
    @inbounds plank.PRO_TN[i]  += ΔT * plank.SP_TN[i]
    @inbounds plank.PRO_TP[i]  += ΔT * plank.SP_TP[i]
    @inbounds plank.PRO_TFe[i] += ΔT * plank.SP_TFe[i]
    @inbounds plank.PRO_RS[i]  += ΔT * plank.SP_RS[i]
    @inbounds plank.PRO_PS[i]  += ΔT * plank.SP_PS[i]
    @inbounds plank.DNA[i]     += ΔT * plank.SDNA[i]
    @inbounds plank.RNA[i]     += ΔT * plank.SRNA[i]

    @inbounds SPRO = plank.SP_RB[i] + plank.SP_MC[i] + plank.SP_MN[i] + plank.SP_TN[i] + 
                     plank.SP_TP[i] + plank.SP_TFe[i]+ plank.SP_RS[i] + plank.SP_PS[i]

    @inbounds plank.CH[i]  -= ΔT * (SPRO + plank.SDNA[i] + plank.SRNA[i])
    @inbounds plank.CH[i]  -= ΔT * (plank.SP_PS[i] * plank.ρChl[i] * 55.0f0 / 893.49f0)

    @inbounds plank.qNH4[i]-= ΔT * (SPRO * p.R_NC_PRO + plank.SDNA[i] * p.R_NC_DNA)
    @inbounds plank.qNH4[i]-= ΔT * (plank.SRNA[i] * p.R_NC_RNA)
    @inbounds plank.qNH4[i]-= ΔT * (plank.SP_PS[i] * plank.ρChl[i] * 4.0f0 / 55.0f0)

    @inbounds plank.PST[i] -= ΔT *(plank.SDNA[i] * p.R_PC_DNA + plank.SRNA[i] * p.R_PC_RNA)
    @inbounds plank.Chl[i] += ΔT * plank.SP_PS[i] * plank.ρChl[i] * 893.49f0 / 55.0f0 # chl unit is mgChl/cell
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
    @inbounds tot_C = total_C_biomass(plank.PRO_RB[i], plank.PRO_MC[i], plank.PRO_MN[i], 
                                      plank.PRO_TN[i], plank.PRO_TP[i], plank.PRO_TFe[i], 
                                      plank.PRO_RS[i], plank.PRO_PS[i], plank.DNA[i], 
                                      plank.RNA[i], plank.CH[i], plank.Chl[i])
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
