##### temperature function for carbon fixation
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

##### calculate ammonium uptake rate (mmolN/individual/second)
@inline function calc_NH4_uptake(NH4, T, CH, qNH4, Bm, pop, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(T, p) * Bm * ac
    return min(VNH4, NH4/ΔT/max(1.0f0,pop)) * p.is_nr
end

##### calculate nitrate uptake rate (mmolN/individual/second)
@inline function calc_NO3_uptake(NO3, T, CH, qNO3, Bm, pop, p, ac, ΔT)
    Qn = qNO3/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNO3max, 1.0f-4, pow = 2.0f0)
    VNO3 = p.VNO3max * regQN * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(T, p) * Bm * ac
    return min(VNO3, NO3/ΔT/max(1.0f0,pop)) * p.is_nr
end

##### calculate phosphate uptake rate (mmolN/individual/second)
@inline function calc_P_uptake(PO4, T, CH, qP, Bm, pop, p, ac, ΔT)
    Qp = (qP + Bm * p.R_PC)/max(1.0f-30, Bm + CH)
    regQP = shape_func_dec(Qp, p.qPmax, 1.0f-4, pow = 2.0f0)
    VPO4 = p.VPO4max * regQP * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(T, p) * Bm * ac
    return min(VPO4, PO4/ΔT/max(1.0f0,pop))
end

##### calculate iron uptake rate (mmolFe/individual/second)
@inline function calc_Fe_uptake(FeT, T, qFe, qFePS, qFeNR, qFeNF, Bm, CH, Sz, pop, p, ac, ΔT)
    Qfe = (qFe + qFePS + qFeNF + qFeNR)/max(1.0f-30, Bm + CH)
    regQFe = shape_func_dec(Qfe, p.qFemax, 1.0f-4, pow = 2.0f0)
    SA = p.SA * Sz^(2.0f0/3.0f0)
    VFe = p.KSAFe * SA * FeT * regQFe * p.Nsuper * tempFunc(T, p) * ac
    return min(VFe, FeT/ΔT/max(1.0f0,pop))
end

@kernel function calc_trs_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VNH4[i] = calc_NH4_uptake(trs.NH4[i], trs.T[i], plank.CH[i], plank.qNH4[i],
                                              plank.Bm[i], trs.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VNO3[i] = calc_NO3_uptake(trs.NO3[i], trs.T[i], plank.CH[i], plank.qNO3[i],
                                              plank.Bm[i], trs.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VPO4[i] = calc_P_uptake(trs.PO4[i], trs.T[i], plank.CH[i], plank.qP[i],
                                            plank.Bm[i], trs.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VFe[i]  = calc_Fe_uptake(trs.FeT[i], trs.T[i], plank.qFe[i], plank.qFePS[i], 
                                             plank.qFeNR[i], plank.qFeNF[i], plank.Bm[i], plank.CH[i], 
                                             plank.Sz[i], trs.pop[i], p, plank.ac[i], ΔT)
end
function calc_trs_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_trs_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### update N, P, Fe quotas - first time
@kernel function update_state_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.qP[i]   += plank.VPO4[i] * ΔT
    @inbounds plank.qNO3[i] += plank.VNO3[i] * ΔT
    @inbounds plank.qNH4[i] += plank.VNH4[i] * ΔT
    @inbounds plank.qFe[i]  += plank.VFe[i]  * ΔT
end
function update_state_1!(plank, ΔT, arch::Architecture)
    kernel! = update_state_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate light reaction (mmolATP/individual/second)
@inline function calc_PS(par, Bm, CH, qFePS, Chl, p)
    αI = par * p.α * Chl / max(1.0f-30, Bm)
    Qfe_ps = qFePS / max(1.0f-30, Bm + CH)
    Ksat = Qfe_ps / max(1.0f-30, Qfe_ps + p.KfePS)
    PS  = p.PCmax * (1.0f0 - exp(-αI)) * Bm * Ksat
    return PS
end
@kernel function calc_PS_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(trs.par[i], plank.Bm[i], plank.CH[i],
                                    plank.qFePS[i], plank.Chl[i], p)
end
function calc_PS!(plank, trs, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate potential maximum respiration (mmolC/individual/second) 
##### and energy production (mmolATP/individual/second)
@inline function calc_respir(CH, Bm, T, p, ac, ΔT)
    RS = CH * p.k_rs * tempFunc(T, p) * ac
    RS = min(RS, CH/ΔT) # double check CH is not over consumed
    ERS = RS * p.e_rs * ac
    return RS, ERS
end
@kernel function calc_respiration_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.RS[i], plank.ERS[i] = calc_respir(plank.CH[i], plank.Bm[i], 
                                                      trs.T[i], p, plank.ac[i], ΔT)
end
function calc_repiration!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_respiration_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate potential carbon fixation rate (mmolC/individual/second)
##### and energy consumption (mmolATP/individual/second)
@inline function calc_CF(CH, Bm, T, p, ac)
    Qc = CH/max(1.0f-30, Bm + CH)
    regQC = shape_func_dec(Qc, p.CHmax, 1.0f-4, pow = 2.0f0)
    CF = p.k_cf * regQC * tempFunc_CF(T, p) * Bm * ac
    ECF = CF * p.e_cf * ac
    return CF, ECF
end
@kernel function calc_carbon_fixation_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.CF[i], plank.ECF[i] = calc_CF(plank.CH[i], plank.Bm[i],
                                                  trs.T[i], p, plank.ac[i])
end
function calc_carbon_fixation!(plank, trs, p, arch::Architecture)
    kernel! = calc_carbon_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate potential nitrate reduction (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@inline function calc_NR(qNO3, qNH4, qFeNR, Bm, CH, T, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    Qfe_NR = qFeNR / max(1.0f-30, Bm + CH)
    Ksat = Qfe_NR / max(1.0f-30, Qfe_NR + p.KfeNR)
    NR = p.k_nr * reg * qNO3 * Ksat * tempFunc(T, p) * ac
    NR = min(NR, qNO3/ΔT) # double check qNO3 are not over consumed
    ENR = NR * p.e_nr * ac
    return NR * p.is_nr, ENR * p.is_nr
end
@kernel function calc_NO3_reduction_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.NR[i], plank.ENR[i] = calc_NR(plank.qNO3[i], plank.qNH4[i],
                                                  plank.qFeNR[i], plank.Bm[i], plank.CH[i],
                                                  trs.T[i], p, plank.ac[i], ΔT)
end
function calc_NO3_reduction!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_NO3_reduction_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate potential nitrogen fixation (mmolN/individual/second)
##### and energy consumption (mmolATP/individual/second)
@inline function calc_NF(qNH4, qFeNF, Bm, CH, T, p, ac)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    Qfe_NF = qFeNF / max(1.0f-30, Bm + CH)
    Ksat = Qfe_NF / max(1.0f-30, Qfe_NF + p.KfeNF)
    NF = p.k_nf * reg * Ksat * tempFunc(T, p) * Bm * ac
    ENF = NF * p.e_nf * ac
    return NF * (p.is_croc + p.is_tric), ENF * (p.is_croc + p.is_tric)
end
@kernel function calc_nitrogen_fixation_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.NF[i], plank.ENF[i] = calc_NF(plank.qNH4[i], plank.qFeNF[i], 
                                                  plank.Bm[i], plank.CH[i],
                                                  trs.T[i], p, plank.ac[i])
end
function calc_nitrogen_fixation!(plank, trs, p, arch::Architecture)
    kernel! = calc_nitrogen_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### energy allocation
@inline function energy_alloc(PS, ERS, ECF, ENF, ENR, p)
    ERSt = ERS; ECFt = ECF; ENFt = ENF; ENRt = ENR;
    exEn = 0.0f0
    if PS ≥ ECF + ENF + ENR
        ERSt = 0.0f0
        exEn = PS - ECF - ENF - ENR
    else #### PS < ECF + ENF + ENR
        if PS + ERS ≥ ECF + ENF + ENR
            ERSt = ECF + ENF + ENR - PS
        else #### PS + ERS < ECF + ENF + ENR
            ECFt = min(PS, ECF)
            ENFt = min(ENF, (PS + ERS - ECFt) * (p.is_croc + p.is_tric))
            ENRt = min(ENR, (PS + ERS - ECFt) * p.is_nr)
        end
    end
    return ERSt, ECFt, ENFt, ENRt, exEn
end
@kernel function energy_allocation_kernel!(plank, p)
    i = @index(Global)
    plank.ERS[i], plank.ECF[i], plank.ENF[i], plank.ENR[i], plank.exEn[i] = 
        energy_alloc(plank.PS[i], plank.ERS[i], plank.ECF[i], 
                     plank.ENF[i], plank.ENR[i], p)
    @inbounds plank.RS[i] = plank.ERS[i] / p.e_rs
    @inbounds plank.CF[i] = plank.ECF[i] / p.e_cf
    @inbounds plank.NF[i] = plank.ENF[i] / p.e_nf
    @inbounds plank.NR[i] = plank.ENR[i] / p.e_nr
end
function energy_allocation!(plank, p, arch::Architecture)
    kernel! = energy_allocation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update energy reserve and CH, and N quotas - second time
@kernel function update_state_2_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.CH[i]   += (plank.CF[i] - plank.RS[i]) * ΔT
    @inbounds plank.qNO3[i] += -plank.NR[i]  * ΔT
    @inbounds plank.qNH4[i] += plank.NR[i]  * ΔT
    @inbounds plank.qNH4[i] += plank.NF[i]  * ΔT
end
function update_state_2!(plank, ΔT, arch::Architecture)
    kernel! = update_state_2_kernel!(device(arch), 256, (size(plank.ac,1)))
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

##### iron allocation
@inline function iron_alloc(par, dpar, qFe, qFePS, qFeNR, qFeNF, tdark, p, ΔT)
    # photosynthesis
    f_ST2PS = p.k_Fe_ST2PS * qFe * isless(0.0f0, dpar)
    f_PS2ST = p.k_Fe_PS2ST * qFePS * (isless(dpar, 0.0f0) + isequal(0.0f0, par))

    # nitrate reduction
    f_ST2NR = p.k_Fe_ST2NR * qFe * isequal(0.0f0, par) * isless(tdark, p.NF_clock) * p.is_nr 
    f_NR2ST = p.k_Fe_NR2ST * qFeNR * (isequal(0.0f0, par) * 
                isless(p.NF_clock, tdark) + isless(0.0f0, par)) * p.is_nr 

    # nitrogen fixation - Crocosphaera watsonii
    f_ST2NF_cr = p.k_Fe_ST2NF * qFe * isequal(0.0f0, par) * isless(tdark, p.NF_clock)
    f_NF2ST_cr = p.k_Fe_NF2ST * qFeNF * (isequal(0.0f0, par) * 
                    isless(p.NF_clock, tdark) + isless(0.0f0, par))
    # nitrogen fixation - Trichodesmium
    f_ST2NF_tr = p.k_Fe_ST2NF * qFe * isless(0.0f0, dpar)
    f_NF2ST_tr = p.k_Fe_NF2ST * qFeNF * (isless(dpar, 0.0f0) + isequal(0.0f0, par))
    # nitrogen fixation
    f_ST2NF = f_ST2NF_cr * p.is_croc + f_ST2NF_tr * p.is_tric
    f_NF2ST = f_NF2ST_cr * p.is_croc + f_NF2ST_tr * p.is_tric

    # photosynthesis
    f_ST2PS = min(f_ST2PS, qFe/ΔT)
    f_PS2ST = min(f_PS2ST, qFePS/ΔT)
    # photosynthesis has the highest priority
    f_ST2NR = min(f_ST2NR, qFe/ΔT - f_ST2PS)
    f_NR2ST = min(f_NR2ST, qFeNR/ΔT)
    # nitrogen fixation
    f_ST2NF = min(f_ST2NF, qFe/ΔT - f_ST2PS)
    f_NF2ST = min(f_NF2ST, qFeNF/ΔT) 

    return f_ST2PS, f_PS2ST, f_ST2NR, f_NR2ST, f_ST2NF, f_NF2ST
end

@kernel function calc_iron_fluxes_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.ST2PS[i], plank.PS2ST[i], plank.ST2NR[i], 
              plank.NR2ST[i], plank.ST2NF[i], plank.NF2ST[i] = 
                    iron_alloc(trs.par[i], trs.dpar[i], plank.qFe[i], 
                    plank.qFePS[i], plank.qFeNR[i], plank.qFeNF[i], plank.tdark[i], p, ΔT)
end
function calc_iron_fluxes!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_iron_fluxes_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
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
    @inbounds plank.qFe[i]  += ΔT * (plank.PS2ST[i] - plank.ST2PS[i] +
                                     plank.NR2ST[i] - plank.ST2NR[i] +
                                     plank.NF2ST[i] - plank.ST2NF[i])
    @inbounds plank.qFePS[i]+= ΔT * (plank.ST2PS[i] - plank.PS2ST[i])
    @inbounds plank.qFeNR[i]+= ΔT * (plank.ST2NR[i] - plank.NR2ST[i])
    @inbounds plank.qFeNF[i]+= ΔT * (plank.ST2NF[i] - plank.NF2ST[i])
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

##### circadian clock after sunset
@kernel function update_tdark_kernel!(plank, trs, ΔT)
    i = @index(Global)
    @inbounds plank.tdark[i] = plank.tdark[i] * (1.0f0 - isless(0.0f0, trs.par[i])) + 
                                ΔT * isequal(0.0f0, trs.par[i])
end
function update_tdark!(plank, trs, ΔT, arch::Architecture)
    kernel! = update_tdark_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, ΔT)
    return nothing
end
