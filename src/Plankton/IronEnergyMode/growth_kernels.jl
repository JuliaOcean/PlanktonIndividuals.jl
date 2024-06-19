##### temperature function for photosynthesis
@inline function tempFunc_PS(temp, p)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0)))*(1.0f0-exp(temp - p.Tmax))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return k/OGT_rate
end

##### temperature function for nutrient uptakes
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0)))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return k/OGT_rate
end

##### intracellular iron allocation to photosynthesis
@inline function iron_alloc_PS(par, p)
    f = shape_func_dec(par, p.Imax, 1.0f-2)
    return f
end

##### calculate light reaction (kJ/individual/second)
@inline function calc_PS(par, Bm, Chl, qFe, p)
    qFe_PS = iron_alloc_PS(par, p) * qFe / max(1.0f-30, Bm)
    αI = par * p.α * Chl / max(1.0f-30, Bm)
    PS  = p.PCmax * (1.0f0 - exp(-αI)) * qFe_PS / max(1.0f-30, qFe_PS + p.KfePS) * Bm
    return PS
end

@kernel function calc_PS_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], plank.Bm[i], plank.Chl[i], plank.qFe[i], p)
end
function calc_PS!(plank, nuts, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### update energy reserve - first time
@kernel function update_state_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i] += plank.PS[i] * ΔT
end
function update_state_1!(plank, ΔT, arch::Architecture)
    kernel! = update_state_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate carbon fixation rate (mmolC/individual/second)
@inline function calc_CF(En, temp, p, ac)
    CF = p.k_cf * En * tempFunc_PS(temp, p) * ac
    ECF = CF * p.e_cf * ac
    return CF, ECF
end

##### calculate respiration (mmolC/individual/second)
@inline function calc_respiration(CH, T, p, ac)
    RS = max(0.0f0, CH) * p.k_rs * tempFunc(T, p) * ac
    ERS = RS * p.e_rs * ac
    return RS, ERS
end

##### carbon fixation and respiration
@kernel function calc_carbon_fluxes_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.CF[i], plank.ECF[i] = calc_CF(plank.En[i], nuts.T[i], p, plank.ac[i])
    @inbounds plank.RS[i], plank.ERS[i] = calc_respiration(plank.CH[i], nuts.T[i], p, plank.ac[i])
end
function calc_carbon_fluxes!(plank, nuts, p, arch::Architecture)
    kernel! = calc_carbon_fluxes_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### calculate ammonium uptake rate (mmolN/individual/second)
@inline function calc_NH4_uptake(NH4, temp, CH, qNH4, Bm, p, ac)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNH4max, 1.0f-4)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(temp, p) * Bm * ac
    return VNH4
end

##### calculate nitrate uptake rate (mmolN/individual/second)
@inline function calc_NO3_uptake(NO3, temp, CH, qNO3, Bm, p, ac)
    Qn = qNO3/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNO3max, 1.0f-4)
    VNO3 = p.VNO3max * regQN * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(temp, p) * Bm * ac
    return VNO3
end

##### calculate phosphate uptake rate (mmolN/individual/second)
@inline function calc_P_uptake(PO4, temp, CH, qP, Bm, p, ac)
    Qp = (qP + Bm * p.R_PC)/max(1.0f-30, Bm + CH)
    regQP = shape_func_dec(Qp, p.qPmax, 1.0f-4)
    VPO4 = p.VPO4max * regQP * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(temp, p) * Bm * ac
    return VPO4
end

##### calculate iron uptake rate (mmolFe/individual/second)
@inline function calc_Fe_uptake(FeT, qFe, Bm, CH, Sz, p, ac)
    Qfe = qFe/max(1.0f-30, Bm +CH)
    regQFe = shape_func_dec(Qfe, p.qFemax, 1.0f-4)
    SA = p.SA * Sz^(2/3)
    VFe = p.KSAFe * SA * FeT * regQFe * ac
    return VFe
end

@kernel function calc_nuts_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.VNH4[i] = calc_NH4_uptake(nuts.NH4[i], nuts.T[i], plank.CH[i], plank.qNH4[i],
                                              plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VNO3[i] = calc_NO3_uptake(nuts.NO3[i], nuts.T[i], plank.CH[i], plank.qNO3[i],
                                              plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VPO4[i] = calc_P_uptake(nuts.PO4[i], nuts.T[i], plank.CH[i], plank.qP[i],
                                            plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VFe[i]  = calc_Fe_uptake(nuts.FeT[i], plank.qFe[i], plank.Bm[i], plank.CH[i], 
                                             plank.Sz[i], p, plank.ac[i])
end
function calc_nuts_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_nuts_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### update energy reserve and CNPFe quotas - second time
@kernel function update_state_2_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i]   += (plank.ERS[i] - plank.ECF[i]) * ΔT
    @inbounds plank.CH[i]   += (plank.CF[i] - plank.RS[i])  * ΔT
    @inbounds plank.RS[i]   += -max(0.0f0, (0.0f0 - plank.CH[i]) / ΔT)
    @inbounds plank.CH[i]   += max(0.0f0, (0.0f0 - plank.CH[i]))
    @inbounds plank.qP[i]   += plank.VPO4[i] * ΔT
    @inbounds plank.qNO3[i] += plank.VNO3[i] * ΔT
    @inbounds plank.qNH4[i] += plank.VNH4[i] * ΔT
    @inbounds plank.qFe[i]  += plank.VFe[i]  * ΔT
end
function update_state_2!(plank, ΔT, arch::Architecture)
    kernel! = update_state_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### nitrate reduction (mmolN/individual/second)
@inline function calc_NO3_reduction(En, qNO3, qNH4, qFe, par, p, ac)
    reg = shape_func_dec(qNH4, p.qNH4max, 1.0f-4)
    qFe_NR = (1.0f0 - iron_alloc_PS(par, p)) * qFe
    NR = p.k_nr * reg * qNO3 * En * qFe_NR / max(1.0f-30, qFe_NR + p.KfeNR) * ac
    ENR = NR * p.e_nr * ac
    return NR, ENR
end

@kernel function calc_NO3_reduc_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.NR[i], plank.ENR[i] = calc_NO3_reduction(plank.En[i], plank.qNO3[i], plank.qNH4[i],
                                                             plank.qFe[i], nuts.par[i], p, plank.ac[i])
end
function calc_NO3_reduc!(plank, nuts, p, arch::Architecture)
    kernel! = calc_NO3_reduc_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### update energy reserve and N quotas - third time
@kernel function update_state_3_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i]   -= plank.ENR[i] * ΔT
    @inbounds plank.qNO3[i] -= plank.NR[i] * ΔT
    @inbounds plank.qNH4[i] += plank.NR[i] * ΔT
end
function update_state_3!(plank, ΔT, arch::Architecture)
    kernel! = update_state_3_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end


##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = p.Chl2N / max(1.0f-30, par[i] * p.α * plank.Chl[i]/plank.Bm[i]) *
                              isless(1.0f-1, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch::Architecture)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, par, p)
    return nothing
end

##### calculate biosynthesis (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.BS[i] = min(plank.CH[i], plank.qNH4[i]/p.R_NC, plank.qP[i]/p.R_PC) * 
                           p.k_mtb * plank.ac[i]
end
function calc_BS!(plank, p, arch::Architecture)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P quotas, biomass, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]   += ΔT * plank.BS[i]
    @inbounds plank.CH[i]   -= ΔT * plank.BS[i]
    @inbounds plank.qNH4[i] -= ΔT * plank.BS[i] * p.R_NC
    @inbounds plank.qP[i]   -= ΔT * plank.BS[i] * p.R_PC
    @inbounds plank.Chl[i]  += ΔT * plank.BS[i] * p.R_NC * plank.ρChl[i]
    @inbounds plank.age[i]  += ΔT / 3600.0f0 * plank.ac[i]
end
function update_biomass!(plank, p, ΔT, arch::Architecture)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### update cell size
##### only calculate functional biomass
@kernel function update_cellsize_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.Sz[i] = plank.Bm[i] / (p.Cquota * p.Nsuper)
end
function update_cellsize!(plank, p, arch::Architecture)
    kernel! = update_cellsize_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end