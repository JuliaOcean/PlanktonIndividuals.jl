##### temperature function
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp+273.15 - p.T⁺))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.T⁺-2)))
    return k/OGT_rate
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, temp, Chl, PRO, DNA, RNA, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * tempFunc(temp, p)
    Bm = functional_C_biomass(PRO, DNA, RNA)
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-10, PCm) * Chl / max(1.0e-10, Bm))) * Bm
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NP_uptake(NH4, NO3, PO4, temp, CH, NST, PST, PRO, DNA, RNA, p, ac)
    Bm = functional_C_biomass(PRO, DNA, RNA)
    C_tot = total_C_biomass(PRO, DNA, RNA, CH)
    N_tot = total_N_biomass(PRO, DNA, RNA, NST, p)
    P_tot = total_P_biomass(DNA, RNA, PST, p)
    regQN = max(0.0, min(1.0, (p.Nqmax - N_tot / max(1.0e-10, C_tot)) / (p.Nqmax - p.Nqmin)))
    regQP = max(0.0, min(1.0, (p.Pqmax - P_tot / max(1.0e-10, C_tot)) / (p.Pqmax - p.Pqmin)))
    VNH4 = p.VNH4max * regQN * NH4/max(1.0e-10, NH4+p.KsatNH4) * tempFunc(temp, p) * Bm * ac
    VNO3 = p.VNO3max * regQN * NO3/max(1.0e-10, NO3+p.KsatNO3) * tempFunc(temp, p) * Bm * ac
    VPO4 = p.VPO4max * regQP * PO4/max(1.0e-10, PO4+p.KsatPO4) * tempFunc(temp, p) * Bm * ac
    return VNH4, VNO3, VPO4
end

@kernel function calc_inorganic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.Chl[i], plank.PRO[i], plank.DNA[i], plank.RNA[i], p) * plank.ac[i]

    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i] = 
                            calc_NP_uptake(nuts.NH4[i], nuts.NO3[i], nut.PO4[i], nuts.T[i],
                                        plank.CH[i], plank.NST[i], plank.PST[i], plank.PRO[i],
                                        plank.DNA[i], plank.RNA[i], p, plank.ac[i])
end
function calc_inorganic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, nuts, p)
    wait(device(arch), event)
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
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 5% of total C acquisition.
@inline function calc_DOC_uptake(DOC, temp, CH, PRO, DNA, RNA, p)
    Bm = functional_C_biomass(PRO, DNA, RNA)
    C_tot = total_C_biomass(PRO, DNA, RNA, CH)
    regQ = max(0.0, min(1.0, (p.Cqmax - CH / max(1.0e-10, C_tot)) / (p.Cqmax - p.Cqmin)))
    VN = p.VDOCmax * regQ * DOC/max(1.0e-10, DOC+p.KsatDOC) * tempFunc(temp, p) * Bm
    return VN
end
@kernel function calc_organic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(nuts.DOC[i], nuts.T[i], plank.CH[i], plank.PRO[i], plank.DNA[i], plank.RNA[i], p) * plank.ac[i]

    @inbounds plank.VDOC[i] = plank.VDOC[i] * isless(0.05, plank.PS[i]/(plank.VDOC[i]+plank.PS[i]))
end
function calc_organic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, nuts, p)
    wait(device(arch), event)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds Bm = functional_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i])
    @inbounds plank.ρChl[i] = plank.PS[i]/max(1.0e-10, plank.PRO[i] + plank.DNA[i] + plank.RNA[i]) * p.Chl2N / 
                             max(1.0e-10, par[i] * p.α * p.Φ * plank.Chl[i]/max(1.0e-10, plank.PRO[i] + plank.DNA[i] + plank.RNA[i])) *
                             isless(1.0e-1, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, par, p)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.resp[i] = p.respir_a * (plank.PRO[i] + plank.DNA[i] + plank.RNA[i])* tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, T, p)
    wait(device(arch), event)
    return nothing
end


##### update C, N, P reserves for the second time of each time step
##### respiration first use Carbohydrate, if it's not enough, use protein then.
@kernel function update_quotas_2_kernel!(plank, ΔT, p)
    i = @index(Global)
    @inbounds plank.CH[i]  = plank.CH[i]  + (plank.VDOC[i] - plank.resp[i]) * ΔT
    @inbounds plank.CH[i]  = plank.CH[i]  + max(0.0, (0.0 - plank.CH[i]))
    @inbounds plank.PRO[i] = plank.PRO[i] - max(0.0, (0.0 - plank.CH[i]))
    @inbounds plank.NST[i] = plank.NST[i] + max(0.0, (0.0 - plank.CH[i])) * p.R_NC_PRO
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, ΔT, p)
    wait(device(arch), event)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, p)
    i  = @index(Global)
    @inbounds limit_PRO = min(plank.CH[i]/(plank.CH[i] + p.k_sat_pro * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_pro * p.Nsuper * p.R_NC_PRO))
    @inbounds limit_DNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_dna * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_dna * p.Nsuper * p.R_NC_DNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_dna * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_rna * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_rna * p.Nsuper * p.R_NC_RNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_rna * p.Nsuper * p.R_PC_RNA))

    @inbounds plank.S_PRO[i] = p.k_pro_a * plank.RNA[i] * limit_PRO
    @inbounds plank.S_DNA[i] = p.k_dna_a * plank.PRO[i] * limit_DNA * isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0)
    @inbounds plank.S_RNA[i] = p.k_rna_a * plank.PRO[i] * limit_RNA
end
function calc_BS!(plank, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, p)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.PRO[i] += ΔT * plank.S_PRO[i]
    @inbounds plank.DNA[i] += ΔT * plank.S_DNA[i]
    @inbounds plank.RNA[i] += ΔT * plank.S_RNA[i]
    @inbounds plank.CH[i]  -= ΔT *(plank.S_PRO[i] + plank.S_DNA[i] + plank.S_RNA[i] + plank.exu[i])
    @inbounds plank.NST[i] -= ΔT *(plank.S_PRO[i] * p.R_NC_PRO + plank.S_DNA[i] * p.R_NC_DNA + plank.S_RNA[i] * p.R_NC_RNA)
    @inbounds plank.PST[i] -= ΔT *(plank.S_DNA[i] * p.R_PC_DNA + plank.S_RNA[i] * p.R_PC_RNA)
    @inbounds plank.Chl[i] += ΔT * plank.S_pro[i] * p.R_NC_PRO * plank.ρChl[i]
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_biomass!(plank, p, ΔT, arch)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, p, ΔT)
    wait(device(arch), event)
    return nothing
end