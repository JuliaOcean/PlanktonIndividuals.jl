##### temperature function for photosynthesis
@inline function tempFunc_PS(T, p)
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
@inline function calc_PS(par, T, Chl, PRO_P, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * tempFunc_PS(T, p)
    PS  = PCm * (1.0f0 - exp(-αI / max(1.0f-30, PCm) * Chl / max(1.0f-30, PRO_P))) * PRO_P
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NP_uptake(NH4, NO3, PO4, T, NST, PST, PRO_P, PRO_B, DNA, RNA, Chl, pop, p, ac, ΔT)
    N_tot = total_N_biomass(PRO_P, PRO_B, DNA, RNA, NST, Chl, p)
    P_tot = total_P_biomass(DNA, RNA, PST, p)
    R_NST = NST / max(1.0f-30, N_tot)
    R_PST = PST / max(1.0f-30, P_tot)
    regQN = shape_func_dec(R_NST, p.NSTmax, 1.0f-4)
    regQP = shape_func_dec(R_PST, p.PSTmax, 1.0f-4)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(T, p) * PRO_B * ac
    VNO3 = p.VNO3max * regQN * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(T, p) * PRO_B * ac
    VPO4 = p.VPO4max * regQP * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(T, p) * PRO_B * ac
    return min(VNH4, NH4/ΔT/max(1.0f0,pop)), 
           min(VNO3, NO3/ΔT/max(1.0f0,pop)), 
           min(VPO4, PO4/ΔT/max(1.0f0,pop))
end

@kernel function calc_inorganic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(trs.par[i], trs.T[i], plank.Chl[i], plank.PRO_P[i], p) * plank.ac[i]

    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i] = 
                            calc_NP_uptake(trs.NH4[i], trs.NO3[i], trs.PO4[i], trs.T[i],
                                        plank.NST[i], plank.PST[i], plank.PRO_P[i], plank.PRO_B[i],
                                        plank.DNA[i], plank.RNA[i], plank.Chl[i], 
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
@inline function calc_DOC_uptake(DOC, T, CH, PRO_P, PRO_B, DNA, RNA, Chl, pop, p, ΔT)
    C_tot = total_C_biomass(PRO_P, PRO_B, DNA, RNA, CH, Chl)
    R_CH = CH / max(1.0f-30, C_tot)
    regQ = shape_func_dec(R_CH, p.CHmax, 1.0f-4)
    VN = p.VDOCmax * regQ * DOC/max(1.0f-30, DOC+p.KsatDOC) * tempFunc(T, p) * PRO_B
    return min(VN, DOC/ΔT/max(1.0f0,pop))
end
@kernel function calc_organic_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(trs.DOC[i], trs.T[i],
                                              plank.CH[i], plank.PRO_P[i], plank.PRO_B[i], plank.DNA[i],
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
    @inbounds plank.ρChl[i] = plank.PS[i] / max(1.0f-30, plank.PRO_P[i]) /
                              max(1.0f-30, par[i] * p.α * p.Φ * plank.Chl[i]/max(1.0f-30, plank.PRO_P[i])) *
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
    @inbounds plank.resp[i] = p.respir * plank.PRO_B[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end


##### update C, N, P reserves for the second time of each time step
@kernel function update_quotas_2_kernel!(plank, ΔT, p)
    i = @index(Global)
    @inbounds plank.CH[i] += (plank.VDOC[i] - plank.resp[i]) * ΔT
    @inbounds plank.CH[i] = max(0.0f0, plank.CH[i])
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, T, p)
    i  = @index(Global)
    @inbounds delta_CH = max(0.0f0, plank.CH[i] - p.CH_min)
    @inbounds limit_PRO_P = min(delta_CH/(delta_CH + p.k_sat_pro * p.Nsuper),
                                plank.NST[i]/(plank.NST[i] + p.k_sat_pro * p.Nsuper * p.R_NC_PRO_P))
    @inbounds limit_PRO_B = min(delta_CH/(delta_CH + p.k_sat_pro * p.Nsuper),
                                plank.NST[i]/(plank.NST[i] + p.k_sat_pro * p.Nsuper * p.R_NC_PRO_B))
    @inbounds limit_DNA = min(delta_CH/(delta_CH + p.k_sat_dna * p.Nsuper),
                              plank.NST[i]/(plank.NST[i] + p.k_sat_dna * p.Nsuper * p.R_NC_DNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_dna * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(delta_CH/(delta_CH + p.k_sat_rna * p.Nsuper),
                              plank.NST[i]/(plank.NST[i] + p.k_sat_rna * p.Nsuper * p.R_NC_RNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_rna * p.Nsuper * p.R_PC_RNA))

    @inbounds S_PRO_P = p.k_pro_P * plank.RNA[i] * limit_PRO_P * tempFunc(T[i], p)
    @inbounds S_PRO_B = p.k_pro_B * plank.RNA[i] * limit_PRO_B * tempFunc(T[i], p)
    @inbounds plank.S_PRO_P[i] = S_PRO_P
    @inbounds plank.S_PRO_B[i] = S_PRO_B
    @inbounds plank.S_DNA[i] = p.k_dna * plank.PRO_B[i] * limit_DNA * tempFunc(T[i], p) *
                               isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0f0)
    @inbounds plank.S_RNA[i] = p.k_rna * plank.PRO_B[i] * limit_RNA * tempFunc(T[i], p)
end
function calc_BS!(plank, T, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds SP = plank.S_PRO_P[i]
    @inbounds SB = plank.S_PRO_B[i]
    @inbounds plank.PRO_P[i] += ΔT * SP
    @inbounds plank.PRO_B[i] += ΔT * SB
    @inbounds plank.DNA[i] += ΔT * plank.S_DNA[i]
    @inbounds plank.RNA[i] += ΔT * plank.S_RNA[i]
    @inbounds plank.CH[i]  -= ΔT *((SP + SB) + plank.S_DNA[i] + plank.S_RNA[i] + 
                                   SP * plank.ρChl[i])
    @inbounds plank.NST[i] -= ΔT *(SP * p.R_NC_PRO_P + SB * p.R_NC_PRO_B + 
                                   plank.S_DNA[i] * p.R_NC_DNA + 
                                   plank.S_RNA[i] * p.R_NC_RNA + 
                                   SP * plank.ρChl[i] * 4.0f0 / 55.0f0)
    @inbounds plank.PST[i] -= ΔT *(plank.S_DNA[i] * p.R_PC_DNA + plank.S_RNA[i] * p.R_PC_RNA)
    @inbounds plank.Chl[i] += ΔT * SP * plank.ρChl[i] * 893.49f0 / 55.0f0 # chl unit is mgChl/cell
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
    @inbounds tot_C = total_C_biomass(plank.PRO_P[i], plank.PRO_B[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i])
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
