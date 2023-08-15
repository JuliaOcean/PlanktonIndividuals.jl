##### temperature function
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp - p.T⁺))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.T⁺+273.15-2)))
    return k/OGT_rate
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, temp, Chl, PRO, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * tempFunc(temp, p)
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-30, PCm) * Chl / max(1.0e-30, PRO))) * PRO
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NP_uptake(NH4, NO3, PO4, temp, NST, PST, PRO, DNA, RNA, Chl, p, ac)
    N_tot = total_N_biomass(PRO, DNA, RNA, NST, Chl, p)
    P_tot = total_P_biomass(DNA, RNA, PST, p)
    R_NST = NST / max(1.0e-30, N_tot)
    R_PST = PST / max(1.0e-30, P_tot)
    fNST = max(0.0, min(1.0, 1.0 - R_NST / p.NSTmax))
    fPST = max(0.0, min(1.0, 1.0 - R_PST / p.PSTmax))
    regQN = fNST^4 / (1.0e-4 + fNST^4)
    regQP = fPST^4 / (1.0e-4 + fPST^4)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0e-30, NH4+p.KsatNH4) * tempFunc(temp, p) * PRO * ac
    VNO3 = p.VNO3max * regQN * NO3/max(1.0e-30, NO3+p.KsatNO3) * tempFunc(temp, p) * PRO * ac
    VPO4 = p.VPO4max * regQP * PO4/max(1.0e-30, PO4+p.KsatPO4) * tempFunc(temp, p) * PRO * ac
    return VNH4, VNO3, VPO4
end

@kernel function calc_inorganic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.Chl[i], plank.PRO[i], p) * plank.ac[i]

    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i] = 
                            calc_NP_uptake(nuts.NH4[i], nuts.NO3[i], nuts.PO4[i], nuts.T[i],
                                        plank.NST[i], plank.PST[i], plank.PRO[i],
                                        plank.DNA[i], plank.RNA[i], plank.Chl[i], p, plank.ac[i])
end
function calc_inorganic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
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
@inline function calc_DOC_uptake(DOC, temp, CH, PRO, DNA, RNA, Chl, p)
    C_tot = total_C_biomass(PRO, DNA, RNA, CH, Chl)
    R_CH = CH / max(1.0e-30, C_tot)
    regQ = max(0.0, min(1.0, (1.0 - R_CH / p.CHmax) * (1.0/(0.01 + 1.0 - R_CH / p.CHmax))))
    # regQ = max(0.0, min(1.0, (p.CHmax - CH / max(1.0e-30, C_tot)) / (p.CHmax - p.CHmin)))
    VN = p.VDOCmax * regQ * DOC/max(1.0e-30, DOC+p.KsatDOC) * tempFunc(temp, p) * PRO
    return VN
end
@kernel function calc_organic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(nuts.DOC[i], nuts.T[i], plank.CH[i], plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.Chl[i], p) * plank.ac[i]

    @inbounds plank.VDOC[i] = plank.VDOC[i] * isless(0.05, plank.PS[i]/(plank.VDOC[i]+plank.PS[i]))
end
function calc_organic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = plank.PS[i]/max(1.0e-30, plank.PRO[i]) / max(1.0e-30, par[i] * p.α * p.Φ * plank.Chl[i]/max(1.0e-30, plank.PRO[i])) *
                             isless(1.0e-1, par[i]) * plank.ac[i]
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
    @inbounds plank.CH[i]  = plank.CH[i]  + max(0.0, (0.0 - plank.CH[i]))
    @inbounds plank.PRO[i] = plank.PRO[i] - max(0.0, (0.0 - plank.CH[i]))
    @inbounds plank.NST[i] = plank.NST[i] + max(0.0, (0.0 - plank.CH[i])) * p.R_NC_PRO
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate protein, DNA, RNA synthesis (mmol C /individual/second)
@kernel function calc_BS_kernel!(plank, T, p)
    i  = @index(Global)
    @inbounds limit_PRO = min(plank.CH[i]/(plank.CH[i] + p.k_sat_pro * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_pro * p.Nsuper * p.R_NC_PRO))
    @inbounds limit_DNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_dna * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_dna * p.Nsuper * p.R_NC_DNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_dna * p.Nsuper * p.R_PC_DNA))
    @inbounds limit_RNA = min(plank.CH[i]/(plank.CH[i] + p.k_sat_rna * p.Nsuper), plank.NST[i]/(plank.NST[i] + p.k_sat_rna * p.Nsuper * p.R_NC_RNA),
                                plank.PST[i]/(plank.PST[i] + p.k_sat_rna * p.Nsuper * p.R_PC_RNA))

    @inbounds plank.S_PRO[i] = p.k_pro * plank.RNA[i] * limit_PRO * tempFunc(T[i], p)
    @inbounds plank.S_DNA[i] = p.k_dna * plank.PRO[i] * limit_DNA * tempFunc(T[i], p) * isless(plank.DNA[i]/(p.C_DNA * p.Nsuper), 2.0)
    @inbounds plank.S_RNA[i] = p.k_rna * plank.PRO[i] * limit_RNA * tempFunc(T[i], p)
end
function calc_BS!(plank, T, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### update C, N, P reserves, protein, DNA, RNA, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds S_Chl = plank.S_PRO[i] * plank.ρChl[i]
    @inbounds plank.PRO[i] += ΔT * plank.S_PRO[i]
    @inbounds plank.DNA[i] += ΔT * plank.S_DNA[i]
    @inbounds plank.RNA[i] += ΔT * plank.S_RNA[i]
    @inbounds plank.CH[i]  -= ΔT *(plank.S_PRO[i] + plank.S_DNA[i] + plank.S_RNA[i] + S_Chl)
    @inbounds plank.NST[i] -= ΔT *(plank.S_PRO[i] * p.R_NC_PRO + plank.S_DNA[i] * p.R_NC_DNA + plank.S_RNA[i] * p.R_NC_RNA + S_Chl*4.0/55.0)
    @inbounds plank.PST[i] -= ΔT *(plank.S_DNA[i] * p.R_PC_DNA + plank.S_RNA[i] * p.R_PC_RNA)
    @inbounds plank.Chl[i] += ΔT * S_Chl * 893.49 / 55.0 # chl unit is mgChl/cell
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_biomass!(plank, p, ΔT, arch)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### calculate exudation of carbon. Nitrogen and phosphorus will not be exuded for now
@kernel function calc_exudation_kernel!(plank, p)
    i = @index(Global)
    @inbounds tot_C = total_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.CH[i], plank.Chl[i])
    @inbounds plank.exu[i] = max(0.0, plank.CH[i] - p.CHmax * tot_C)
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
