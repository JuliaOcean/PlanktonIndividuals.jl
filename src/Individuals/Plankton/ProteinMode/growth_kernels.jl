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
@inline function calc_PS(par, PRO_P, Chl, p)
    αI = par * p.α * Chl / max(1.0f-30, PRO_P)
    PS = p.PCmax * (1.0f0 - exp(-αI)) * PRO_P
    return PS
end

##### calculate photosynthesis rate kernel
@kernel function calc_PS_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(trs.par[i], plank.PRO_P[i], plank.Chl[i], p) * plank.ac[i]
end
function calc_PS!(plank, trs, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@inline function calc_respir(CH, PRO_RS, T, p, ac, ΔT)
    RS = p.KcatRS * PRO_RS * tempFunc(T, p) * ac
    RS = min(RS, CH/ΔT)
    ERS = RS * p.e_rs * ac
    return RS, ERS
end
@kernel function calc_respir_kernel!(plank, T, p, ΔT)
    i = @index(Global)
    @inbounds plank.RS[i], plank.ERS[i] = calc_respir(plank.CH[i], plank.PRO_RS[i], T[i], p, plank.ac[i], ΔT)

end
function calc_respiration!(plank, T, p, ΔT, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p, ΔT)
    return nothing
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NPFe_uptake(NH4, NO3, PO4, DFe, T, PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, qNH3, qNO3, PST, qFe, DNA, RNA, CH, Chl, pop, p, ac, ΔT)
    C_tot = total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl)
    N_tot = total_N_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, qNO3, qNH3, Chl, p)
    P_tot = total_P_biomass(DNA, RNA, PST, p)
    R_qNH3 = qNH3 / max(1.0f-30, N_tot)
    R_qNO3 = qNO3 / max(1.0f-30, N_tot)
    R_PST  = PST  / max(1.0f-30, P_tot)
    R_qFe  = qFe  / max(1.0f-30, C_tot)
    regQNH3 = shape_func_dec(R_qNH3, p.qNH4max, 1.0f-4)
    regQNO3 = shape_func_dec(R_qNO3, p.qNO3max, 1.0f-4)
    regQP   = shape_func_dec(R_PST,  p.PSTmax,  1.0f-4)
    regQFe  = shape_func_dec(R_qFe,  p.qFemax,  1.0f-4)
    VNH4 = min( p.KcatTNH4 * PRO_Tn  * NH4 / max(1.0f-30, NH4 + p.KsatNH4) * tempFunc(T, p) * regQNH3  * ac , NH4/ΔT/max(1.0f0,pop))  #mmolN/individual/second
    VNO3 = min( p.KcatTNO3 * PRO_Tn  * NO3 / max(1.0f-30, NO3 + p.KsatNO3) * tempFunc(T, p) * regQNO3  * ac , NO3/ΔT/max(1.0f0,pop))  #mmolN/individual/second    
    VPO4 = min( p.KcatTPO4 * PRO_Tp  * PO4 / max(1.0f-30, PO4 + p.KsatPO4) * tempFunc(T, p) * regQP    * ac , PO4/ΔT/max(1.0f0,pop))  #mmolP/individual/second
    VFe  = min( p.KcatTFe  * PRO_Tfe * DFe / max(1.0f-30, DFe + p.KsatFe)  * tempFunc(T, p) * regQFe   * ac , DFe/ΔT/max(1.0f0,pop))  #mmolFe/individual/second
    
    EVNO3 = VNO3 * p.e_tno3 
    EVPO4 = VPO4 * p.e_tpo4 
    EVFe  = VFe  * p.e_tfe  

    return VNH4, VNO3, VPO4, VFe, EVNO3, EVPO4, EVFe
    
end

@kernel function calc_inorganic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i], plank.VFe[i], plank.EVNO3[i], plank.EVPO4[i], plank.EVFe[i] = 
                            calc_NPFe_uptake(trs.NH4[i], trs.NO3[i], trs.PO4[i], trs.DFe[i], trs.T[i], plank.PRO_R[i], plank.PRO_Mc[i], 
                                             plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i],
                                             plank.qNH3[i], plank.qNO3[i], plank.PST[i], plank.qFe[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i],
                                             trs.pop[i], p, plank.ac[i], ΔT)
end
function calc_inorganic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate energy consumption of nutrient uptakes (mmolATP/individual/second)
@inline function calc_NPFe_energy(PS, ERS, EVNO3, EVPO4, EVFe, p)
    ERS_avail = max(ERS-p.e_min, 0.0f0)
    EPS_avail = max(PS - max(p.e_min - ERS, 0.0f0), 0.0f0)
    Esupply = ERS_avail + EPS_avail
    Edemand = EVNO3 + EVPO4 + EVFe
    if Edemand <= 0.0f0
        return 0.0f0, 0.0f0, 0.0f0, ERS_avail, EPS_avail
    else
        Eused = min(Esupply, Edemand)
        scale = Eused / Edemand
        tEVNO3 = EVNO3 * scale
        tEVPO4 = EVPO4 * scale
        tEVFe  = EVFe  * scale
        exEN_RS = max(ERS_avail - Eused, 0.0f0)
        exEN_PS = max(EPS_avail - max(Eused - ERS_avail, 0.0f0), 0.0f0)
    end
    return tEVNO3, tEVPO4, tEVFe, exEN_RS, exEN_PS
end

@kernel function calc_NPFe_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.EVNO3[i], plank.EVPO4[i], plank.EVFe[i], plank.exEN_RS[i], plank.exEN_PS[i] = 
        calc_NPFe_energy(plank.PS[i], plank.ERS[i], plank.EVNO3[i], plank.EVPO4[i], plank.EVFe[i], p)
       
    @inbounds plank.VNO3[i] = plank.EVNO3[i] / p.e_tno3 
    @inbounds plank.VPO4[i] = plank.EVPO4[i] / p.e_tpo4
    @inbounds plank.VFe[i]  = plank.EVFe[i]  / p.e_tfe
end

function calc_NPFe_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_NPFe_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P reserves for the first time of each time step
@kernel function update_quotas_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.PST[i]  += plank.VPO4[i] * ΔT
    @inbounds plank.qNH3[i] += plank.VNH4[i] * ΔT
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
@inline function calc_CF(CH, Chl, DNA, RNA, PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, T, p, ac)
    C_tot = total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl)
    Qc    = CH / max(1.0f-30, C_tot)
    regQC = shape_func_dec(Qc, p.CHmax, 1.0f-4)
    CF  = p.KcatCF * PRO_Mc * tempFunc_CF(T, p) * regQC * ac
    ECF = CF * p.e_cf * ac
    return CF, ECF
end
@kernel function calc_carbon_fixation_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.CF[i], plank.ECF[i] = calc_CF(plank.CH[i], plank.Chl[i], plank.DNA[i], plank.RNA[i], plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], 
                                                  plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i], T[i], p, plank.ac[i])
end
function calc_carbon_fixation!(plank, T, p, arch)
    kernel! = calc_carbon_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate potential nitrate reduction (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@inline function calc_NR(qNO3, qNH4, PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl, T, p, ac, ΔT)
    C_tot = total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl)
    Qn = qNH4/max(1.0f-30, C_tot)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    Ksat = qNO3 / max(1.0f-30, (qNO3 + p.KsatNR))

    NR = p.KcatNR * PRO_Mn * tempFunc(T, p) * reg * Ksat * ac
    NR = min(NR, qNO3/ΔT) # double check qNO3 are not over consumed
    ENR = NR * p.e_nr * ac
    return NR * p.is_nr, ENR * p.is_nr
end
@kernel function calc_NO3_reduction_kernel!(plank, T, p, ΔT)
    i = @index(Global)
    @inbounds plank.NR[i], plank.ENR[i] = calc_NR(plank.qNO3[i], plank.qNH4[i], plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i], T[i], p, plank.ac[i], ΔT)
end
function calc_NO3_reduction!(plank, T, p, ΔT, arch::Architecture)
    kernel! = calc_NO3_reduction_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p, ΔT)
    return nothing
end

##### calculate potential nitrogen fixation (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@inline function calc_NF(qNH4, PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl, T, p, ac)
    C_tot = total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl)
    Qn = qNH4/max(1.0f-30, C_tot)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    NF = p.KcatNF * PRO_Mn * tempFunc(T, p) * reg * ac
    ENF = NF * p.e_nf * ac
    return NF * (p.is_croc + p.is_tric), ENF * (p.is_croc + p.is_tric)
end
@kernel function calc_nitrogen_fixation_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.NF[i], plank.ENF[i] = calc_NF(plank.qNH4[i], plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i], T[i], p, plank.ac[i])
end
function calc_nitrogen_fixation!(plank, T, p, arch::Architecture)
    kernel! = calc_nitrogen_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate energy consumption of carbon and nitrigen metabolism (mmolATP/individual/second)
@inline function calc_CN_energy(ECF, ENR, ENF, exEN_PS)
    Esupply = exEN_PS 
    Edemand = ECF + ENR + ENF
    if Edemand <= 0.0f0
        return 0.0f0, 0.0f0, 0.0f0, exEN_PS
    else
        Eused = min(Esupply, Edemand)
        scale = Eused / Edemand
        tECF = ECF * scale
        tENR = ENR * scale
        tENF = ENF * scale
        exEN_PS = max(exEN_PS - Eused, 0.0f0)
    end
        
    return tECF, tENR, tENF, exEN_PS
end
@kernel function calc_CN_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.ECF[i], plank.ENR[i], plank.ENF[i], plank.exEN_PS[i] = 
        calc_CN_energy(plank.ECF[i], plank.ENR[i], plank.ENF[i], plank.exEN_PS[i])
       
    @inbounds plank.CF[i] = plank.ECF[i] / p.e_cf
    @inbounds plank.NR[i] = plank.ENR[i] / p.e_nr
    @inbounds plank.NF[i] = plank.ENF[i] / p.e_nf 
end
function calc_CN_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_CN_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 5% of total C acquisition.
@inline function calc_DOC_uptake(DOC, T, CH, PRO_R, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_Mc, PRO_Mn, PRO_RS, PRO_P, DNA, RNA, Chl, qNO3, qNH3, PST, pop, p, ΔT)
    C_tot = total_C_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, CH, Chl)
    N_tot = total_N_biomass(PRO_R, PRO_Mc, PRO_Mn, PRO_Tn, PRO_Tp, PRO_Tfe, PRO_RS, PRO_P, DNA, RNA, qNO3, qNH3, Chl, p)
    P_tot = total_P_biomass(DNA, RNA, PST, p)
    R_CH = CH / max(1.0f-30, C_tot)
    regQ = shape_func_dec(R_CH, p.CHmax, 1.0f-4)
    VN = p.VDOCmax * regQ * DOC/max(1.0f-30, DOC + p.KsatDOC) * tempFunc(T, p) * min(C_tot, N_tot * 106.0f0 ./ 16.0f0, P_tot * 106.0f0) 
    return min(VN, DOC/ΔT/max(1.0f0,pop))
end

@kernel function calc_organic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(trs.DOC[i], trs.T[i],
                                              plank.CH[i], plank.PRO_R[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i],
                                              plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_RS[i], plank.PRO_P[i], plank.DNA[i], plank.RNA[i], plank.Chl[i],
                                              plank.qNO3[i], plank.qNH3[i], plank.PST[i], trs.pop[i], p, ΔT) * plank.ac[i]

    @inbounds plank.VDOC[i] = plank.VDOC[i] * isless(0.05f0, plank.CF[i]/(plank.VDOC[i]+plank.CF[i]))
end
function calc_organic_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = p.Chl2N / max(1.0f-30, par[i] * p.α * plank.Chl[i] / plank.PRO_P[i]) *
                              isless(1.0f-1, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, par, p)
    return nothing
end

##### calculate protein and chla degration under N stravation (mmolC/individual/second)
@inline function calc_degradation(PRO_R, PRO_P,PRO_Mc, Chl, qNH3, qNO3, p , T, ac)
    NST = qNH3 + qNO3
    reg_Pr  = shape_func_dec(NST, p.NSTmax, 1.0f-4, pow = 1.0f0)
    reg_Pp  = shape_func_dec(NST, p.NSTmax, 1.0f-4, pow = 4.0f0)
    reg_Pmc = shape_func_dec(NST, p.NSTmax, 1.0f-4, pow = 2.0f0)
    reg_chl = shape_func_dec(NST, p.NSTmax, 1.0f-4, pow = 4.0f0)

    D_Pr  = p.k_degr   * (PRO_R  - p.PRO_Rmin)   * reg_Pr  * tempFunc(T, p) * ac
    D_Pp  = p.k_degp   * (PRO_P  - p.PRO_Pmin)   * reg_Pp  * tempFunc(T, p) * ac
    D_Pmc = p.k_degmc  * (PRO_Mc - p.PRO_Mcmin)  * reg_Pmc * tempFunc(T, p) * ac
    D_chl = p.k_degchl * (Chl - p.Chlmin)        * reg_chl * tempFunc(T, p) * ac
    return D_Pr, D_Pp, D_Pmc, D_chl
end
@kernel function calc_degradation_kernel!(plank, p, trs)
    i = @index(Global)
    @inbounds plank.D_Pr[i], plank.D_Pp[i], plank.D_Pmc[i], plank.D_chl[i] = 
        calc_degradation(plank.PRO_R[i], plank.PRO_P[i], plank.PRO_Mc[i], plank.Chl[i], plank.qNH3[i], plank.qNO3[i], p, trs.T[i], plank.ac[i])
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
    @inbounds plank.CH[i]     = plank.CH[i]   + (plank.CF[i] + plank.VDOC[i] ) * ΔT +
                               (plank.D_Pr[i] + plank.D_Pp[i] + plank.D_Pmc[i]) * ΔT + plank.D_chl[i] / 893.49f0 * 55.0f0 * ΔT
    @inbounds plank.qNH3[i]   = plank.qNH3[i] + (plank.NR[i] + plank.NF[i]) * ΔT + 
                               (plank.D_Pr[i] + plank.D_Pp[i] + plank.D_Pmc[i]) * p.R_NC_PRO * ΔT + plank.D_chl[i] / 893.49f0 * 4.0f0 * ΔT
    @inbounds plank.qNO3[i]   = plank.qNO3[i]   - plank.NR[i]    * ΔT
    @inbounds plank.PRO_R[i]  = plank.PRO_R[i]  - plank.D_Pr[i]  * ΔT
    @inbounds plank.PRO_P[i]  = plank.PRO_P[i]  - plank.D_Pp[i]  * ΔT
    @inbounds plank.PRO_Mc[i] = plank.PRO_Mc[i] - plank.D_Pmc[i] * ΔT
    @inbounds plank.Chl[i]    = plank.Chl[i]    - plank.D_chl[i] * ΔT
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, trs, p)
    i  = @index(Global)
    @inbounds NST = plank.qNH3[i] + plank.qNO3[i]
    @inbounds limit_PRO = min(plank.CH[i]/(plank.CH[i] + p.k_sat_pro * p.Nsuper),
                              NST/(NST + p.k_sat_pro * p.Nsuper * p.R_NC_PRO))
    @inbounds limit_DNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_dna * p.Nsuper),
                              NST/(NST + p.k_sat_dna * p.Nsuper * p.R_NC_DNA),
                              plank.PST[i]/(plank.PST[i] + p.k_sat_dna * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_rna * p.Nsuper),
                              NST/(NST + p.k_sat_rna * p.Nsuper * p.R_NC_RNA),
                              plank.PST[i]/(plank.PST[i] + p.k_sat_rna * p.Nsuper * p.R_PC_RNA))
    
    @inbounds C_tot = total_C_biomass(plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i])
    @inbounds QFe  = plank.qFe[i] / max(1.0f-30, C_tot)
    @inbounds Ksat_Fe = QFe / max(1.0f-30, (QFe + p.KFe_em))
    @inbounds αI = trs.par[i] * p.α * plank.Chl[i] / max(1.0f-30, plank.PRO_P[i])
    @inbounds RegLight = αI / (αI + 0.1f0)
    @inbounds lim_TN =  shape_func_dec(trs.NO3[i] + trs.NH4[i], p.TN_max, 1.0f-4)
    @inbounds lim_TP =  shape_func_dec(trs.PO4[i], p.TP_max, 1.0f-4)
    @inbounds lim_TFe = shape_func_dec(trs.DFe[i], p.TFe_max, 1.0f-4)
     
    @inbounds plank.S_Pr[i]  = plank.PRO_R[i]  * p.r_max / p.n_r    * limit_PRO * tempFunc(trs.T[i], p) * limit_RNA * plank.ac[i]
    @inbounds plank.S_Pmc[i] = plank.PRO_R[i]  * p.r_max / p.n_mC   * limit_PRO * tempFunc(trs.T[i], p) * plank.ac[i]
    @inbounds plank.S_Pmn[i] = plank.PRO_R[i]  * p.r_max / p.n_mN   * limit_PRO * tempFunc(trs.T[i], p) * Ksat_Fe   * plank.ac[i]
    @inbounds plank.S_Ptn[i] = plank.PRO_R[i]  * p.r_max / p.n_tN   * limit_PRO * tempFunc(trs.T[i], p) * lim_TN    * plank.ac[i]
    @inbounds plank.S_Ptp[i] = plank.PRO_R[i]  * p.r_max / p.n_tPO4 * limit_PRO * tempFunc(trs.T[i], p) * lim_TP    * plank.ac[i]
    @inbounds plank.S_Ptfe[i]= plank.PRO_R[i]  * p.r_max / p.n_tFe  * limit_PRO * tempFunc(trs.T[i], p) * lim_TFe   * plank.ac[i]
    @inbounds plank.S_Prs[i] = plank.PRO_R[i]  * p.r_max / p.n_rs   * limit_PRO * tempFunc(trs.T[i], p) * Ksat_Fe   * plank.ac[i]
    @inbounds plank.S_Pp[i] = plank.PRO_R[i] * p.r_max / p.n_p * limit_PRO * tempFunc(trs.T[i], p) * RegLight * Ksat_Fe * plank.ac[i]
    @inbounds plank.S_DNA[i] = p.k_dna  * limit_DNA * tempFunc(trs.T[i], p) * isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0f0) * plank.ac[i]
    @inbounds plank.S_RNA[i] = plank.S_Pr[i] * p.R_C_RNAPr * plank.ac[i] 

    @inbounds plank.ESPr[i]  = plank.S_Pr[i] * p.e_sp
    @inbounds plank.ESPmc[i] = plank.S_Pmc[i] * p.e_sp
    @inbounds plank.ESPmn[i] = plank.S_Pmn[i] * p.e_sp
    @inbounds plank.ESPtn[i] = plank.S_Ptn[i] * p.e_sp
    @inbounds plank.ESPtp[i] = plank.S_Ptp[i] * p.e_sp
    @inbounds plank.ESPtfe[i]= plank.S_Ptfe[i] * p.e_sp
    @inbounds plank.ESPrs[i] = plank.S_Prs[i] * p.e_sp
    @inbounds plank.ESPp[i] = plank.S_Pp[i] * p.e_sp
    @inbounds plank.EDNA[i] = plank.S_DNA[i] * p.e_dna
    @inbounds plank.ERNA[i] = plank.S_RNA[i] * p.e_rna  
end
function calc_BS!(plank, trs, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate energy consumption of protein, DNA, RNA synthesis (mmolATP/individual/second)
@inline function calc_PRO_DNA_RNA_energy(ESPr,ESPmc, ESPmn, ESPtn, ESPtp, ESPtfe, ESPrs, ESPp, EDNA, ERNA, exEN_PS, exEN_RS)
    Esupply = exEN_PS + exEN_RS
    PS_avail = exEN_PS
    RS_avail = exEN_RS
    Edemand = ESPr + ESPmc + ESPmn + ESPtn + ESPtp + ESPtfe + ESPrs + ESPp + EDNA + ERNA
    if Edemand <= 0.0f0
        return 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, exEN_PS, exEN_RS
    else
        Eused = min(Esupply, Edemand)
        scale = Eused / Edemand
        tESPr  = ESPr  * scale
        tESPmc = ESPmc * scale
        tESPmn = ESPmn * scale
        tESPtn = ESPtn * scale
        tESPtp = ESPtp * scale
        tESPtfe = ESPtfe * scale
        tESPrs = ESPrs * scale
        tESPp = ESPp * scale
        tEDNA = EDNA * scale
        tERNA = ERNA * scale
        exEN_PS = max(PS_avail - Eused, 0.0f0)
        exEN_RS = max(RS_avail - max(Eused - PS_avail, 0.0f0), 0.0f0)
    end
        
    return tESPr, tESPmc, tESPmn, tESPtn, tESPtp, tESPtfe, tESPrs, tESPp, tEDNA, tERNA, exEN_PS, exEN_RS
end
@kernel function calc_PRO_DNA_RNA_energy_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.ESPr[i], plank.ESPmc[i], plank.ESPmn[i], plank.ESPtn[i], plank.ESPtp[i], plank.ESPtfe[i], plank.ESPrs[i], plank.ESPp[i], plank.EDNA[i], plank.ERNA[i], plank.exEN_PS[i], plank.exEN_RS[i] = 
        calc_PRO_DNA_RNA_energy(plank.ESPr[i], plank.ESPmc[i], plank.ESPmn[i], plank.ESPtn[i], plank.ESPtp[i], plank.ESPtfe[i], plank.ESPrs[i], plank.ESPp[i], plank.EDNA[i], plank.ERNA[i], plank.exEN_PS[i], plank.exEN_RS[i])
        plank.S_Pr[i]  = plank.ESPr[i]  / p.e_sp
        plank.S_Pmc[i] = plank.ESPmc[i]  / p.e_sp
        plank.S_Pmn[i] = plank.ESPmn[i]  / p.e_sp
        plank.S_Ptn[i] = plank.ESPtn[i]  / p.e_sp
        plank.S_Ptp[i] = plank.ESPtp[i]  / p.e_sp
        plank.S_Ptfe[i]= plank.ESPtfe[i] / p.e_sp
        plank.S_Prs[i] = plank.ESPrs[i]  / p.e_sp
        plank.S_Pp[i] = plank.ESPp[i] / p.e_sp
        plank.S_DNA[i] = plank.EDNA[i]   / p.e_dna
        plank.S_RNA[i] = plank.ERNA[i]   / p.e_rna
end
function calc_PRO_DNA_RNA_energy_alloc!(plank, p, arch::Architecture)
    kernel! = calc_PRO_DNA_RNA_energy_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.PRO_R[i]   += ΔT * plank.S_Pr[i]
    @inbounds plank.PRO_Mc[i]  += ΔT * plank.S_Pmc[i]
    @inbounds plank.PRO_Mn[i]  += ΔT * plank.S_Pmn[i]
    @inbounds plank.PRO_Tn[i]  += ΔT * plank.S_Ptn[i]
    @inbounds plank.PRO_Tp[i]  += ΔT * plank.S_Ptp[i]
    @inbounds plank.PRO_Tfe[i] += ΔT * plank.S_Ptfe[i]
    @inbounds plank.PRO_RS[i]  += ΔT * plank.S_Prs[i]
    @inbounds plank.PRO_P[i]   += ΔT * plank.S_Pp[i]
    @inbounds plank.DNA[i] += ΔT * plank.S_DNA[i]
    @inbounds plank.RNA[i] += ΔT * plank.S_RNA[i]

    @inbounds S_PRO = plank.S_Pr[i] + plank.S_Pmc[i] + plank.S_Pmn[i] + plank.S_Ptn[i] + plank.S_Ptp[i] + plank.S_Ptfe[i] + plank.S_Prs[i] + plank.S_Pp[i]
    @inbounds plank.CH[i]  -= ΔT * (S_PRO + plank.S_DNA[i] + plank.S_RNA[i] +
                                    plank.S_Pp[i] * plank.ρChl[i] * 55.0f0 / 893.49f0)
    @inbounds plank.qNH3[i]-= ΔT *(S_PRO * p.R_NC_PRO + plank.S_DNA[i] * p.R_NC_DNA + 
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
    @inbounds tot_C = total_C_biomass(plank.PRO_R[i], plank.PRO_Mc[i], plank.PRO_Mn[i], plank.PRO_Tn[i], plank.PRO_Tp[i], plank.PRO_Tfe[i], plank.PRO_RS[i], plank.PRO_P[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i])
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
