##### temperature function for photosynthesis
@inline function tempFunc_PS(temp, p)
    x = temp - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 4.0f-2)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0))) * regT
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### temperature function for nutrient uptakes
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0)))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### calculate light reaction (kJ/individual/second)
@inline function calc_PS(par, Bm, CH, Chl, qFePS, p)
    Qfe_ps = qFePS / max(1.0f-30, Bm + CH)
    αI = par * p.α * Chl / max(1.0f-30, Bm)
    PS  = p.PCmax * (1.0f0 - exp(-αI)) * Qfe_ps / max(1.0f-30, Qfe_ps + p.KfePS) * Bm
    return PS
end

@kernel function calc_PS_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], plank.Bm[i], plank.CH[i], plank.Chl[i], plank.qFe[i], p)
end
function calc_PS!(plank, nuts, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@inline function calc_respir(CH, En, T, p, ac, ΔT)
    reg = shape_func_dec_alt(En, p.Enmax, 1.0f-2)
    RS = max(0.0f0, CH) * reg * p.k_rs * tempFunc(T, p) * ac
    RS = min(RS, CH/ΔT) # double check CH is not over consumed
    ERS = RS * p.e_rs * ac
    return RS, ERS
end

@kernel function calc_respiration_kernel!(plank, nuts, p, ΔT)
    i = @index(Global)
    @inbounds plank.RS[i], plank.ERS[i] = calc_respir(plank.CH[i], plank.En[i], 
                                                      nuts.T[i], p, plank.ac[i], ΔT)
end
function calc_repiration!(plank, nuts, p, ΔT, arch::Architecture)
    kernel! = calc_respiration_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p, ΔT)
    return nothing
end

##### update energy reserve - first time
@kernel function update_state_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i] += (plank.PS[i] + plank.ERS[i]) * ΔT
    @inbounds plank.CH[i] += -plank.RS[i] * ΔT
end
function update_state_1!(plank, ΔT, arch::Architecture)
    kernel! = update_state_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate carbon fixation rate (mmolC/individual/second)
@inline function calc_CF(En, CH, Bm, temp, p, ac, ΔT)
    Qc = CH/max(1.0f-30, Bm + CH)
    regQC = shape_func_dec(Qc, p.CHmax, 1.0f-4)
    regEn = shape_func_inc_alt(En, p.Enmax, 1.0f-4)
    CF = p.k_cf * regQC * regEn * tempFunc_PS(temp, p) * Bm * ac
    CF = min(CF, En/p.e_cf/ΔT)
    ECF = CF * p.e_cf * ac
    return CF, ECF
end

@kernel function calc_carbon_fixation_kernel!(plank, nuts, p, ΔT)
    i = @index(Global)
    @inbounds plank.CF[i], plank.ECF[i] = calc_CF(plank.En[i], plank.CH[i], plank.Bm[i],
                                                  nuts.T[i], p, plank.ac[i], ΔT)
end
function calc_carbon_fixation!(plank, nuts, p, ΔT, arch::Architecture)
    kernel! = calc_carbon_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p, ΔT)
    return nothing
end

##### calculate ammonium uptake rate (mmolN/individual/second)
@inline function calc_NH4_uptake(NH4, temp, CH, qNH4, Bm, pop, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNH4max, 1.0f-4)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(temp, p) * Bm * ac
    return min(VNH4, NH4/ΔT/max(1.0f0,pop)) * p.is_nr
end

##### calculate nitrate uptake rate (mmolN/individual/second)
@inline function calc_NO3_uptake(NO3, temp, CH, qNO3, Bm, pop, p, ac, ΔT)
    Qn = qNO3/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNO3max, 1.0f-4)
    VNO3 = p.VNO3max * regQN * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(temp, p) * Bm * ac
    return min(VNO3, NO3/ΔT/max(1.0f0,pop)) * p.is_nr
end

##### calculate phosphate uptake rate (mmolN/individual/second)
@inline function calc_P_uptake(PO4, temp, CH, qP, Bm, pop, p, ac, ΔT)
    Qp = (qP + Bm * p.R_PC)/max(1.0f-30, Bm + CH)
    regQP = shape_func_dec(Qp, p.qPmax, 1.0f-4)
    VPO4 = p.VPO4max * regQP * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(temp, p) * Bm * ac
    return min(VPO4, PO4/ΔT/max(1.0f0,pop))
end

##### calculate iron uptake rate (mmolFe/individual/second)
@inline function calc_Fe_uptake(FeT, qFe, Bm, CH, Sz, pop, p, ac, ΔT)
    Qfe = qFe/max(1.0f-30, Bm +CH)
    regQFe = shape_func_dec(Qfe, p.qFemax, 1.0f-4)
    SA = p.SA * Sz^(2/3)
    VFe = p.KSAFe * SA * FeT * regQFe * ac
    return min(VFe, FeT/ΔT/max(1.0f0,pop))
end

@kernel function calc_nuts_uptake_kernel!(plank, nuts, p, ΔT)
    i = @index(Global)
    @inbounds plank.VNH4[i] = calc_NH4_uptake(nuts.NH4[i], nuts.T[i], plank.CH[i], plank.qNH4[i],
                                              plank.Bm[i], nuts.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VNO3[i] = calc_NO3_uptake(nuts.NO3[i], nuts.T[i], plank.CH[i], plank.qNO3[i],
                                              plank.Bm[i], nuts.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VPO4[i] = calc_P_uptake(nuts.PO4[i], nuts.T[i], plank.CH[i], plank.qP[i],
                                            plank.Bm[i], nuts.pop[i], p, plank.ac[i], ΔT)
    @inbounds plank.VFe[i]  = calc_Fe_uptake(nuts.FeT[i], plank.qFe[i], plank.Bm[i], plank.CH[i], 
                                             plank.Sz[i], nuts.pop[i], p, plank.ac[i], ΔT)
end
function calc_nuts_uptake!(plank, nuts, p, ΔT, arch::Architecture)
    kernel! = calc_nuts_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p, ΔT)
    return nothing
end

##### update energy reserve and CNPFe quotas - second time
##### double check energy reserve is non-negative
##### use all the energy if it's not enough
@kernel function update_state_2_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i]   += -plank.ECF[i] * ΔT
    @inbounds plank.CH[i]   += plank.CF[i]   * ΔT
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
@inline function calc_NO3_reduction(En, qNO3, qNH4, qFeNR, Bm, CH, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4)
    Qfe_NR = qFeNR / max(1.0f-30, Bm + CH)
    regEn = shape_func_inc_alt(En, p.Enmax * p.Nsuper, 1.0f-4)
    NR = p.k_nr * reg * regEn * qNO3 * Qfe_NR / max(1.0f-30, Qfe_NR + p.KfeNR) * ac
    NR = min(NR, qNO3/ΔT) # double check qNO3 are not over consumed
    ENR = NR * p.e_nr * ac
    return NR * p.is_nr, ENR * p.is_nr
end

@kernel function calc_NO3_reduc_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.NR[i], plank.ENR[i] = calc_NO3_reduction(plank.En[i], plank.qNO3[i], plank.qNH4[i],
                                                             plank.qFeNR[i], plank.Bm[i], plank.CH[i],
                                                             p, plank.ac[i], ΔT)
end
function calc_NO3_reduc!(plank, p, ΔT, arch::Architecture)
    kernel! = calc_NO3_reduc_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### nitrogen fixation (mmolN/individual/second)
@inline function calc_Nfixation(En, qNH4, qFeNF, Bm, CH, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4)
    Qfe_NF = qFeNF / max(1.0f-30, Bm + CH)
    regEn = shape_func_inc_alt(En, p.Enmax * p.Nsuper, 1.0f-4)
    NF = p.k_nf * reg * regEn * Qfe_NF / max(1.0f-30, Qfe_NF + p.KfeNF) * Bm * ac
    NF = min(NF, En/p.e_nf/ΔT)
    ENF = NF * p.e_nf * ac
    return NF * (p.is_croc + p.is_tric), ENF * (p.is_croc + p.is_tric)
end

@kernel function calc_Nfix_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.NF[i], plank.ENF[i] = calc_Nfixation(plank.En[i], plank.qNH4[i],
                                                         plank.qFeNF[i], plank.Bm[i], plank.CH[i],
                                                         p, plank.ac[i], ΔT)
end
function calc_Nfix!(plank, p, ΔT, arch::Architecture)
    kernel! = calc_Nfix_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### update energy reserve and N quotas - third time
@kernel function update_state_3_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.En[i]   -= plank.ENR[i] * ΔT
    @inbounds plank.En[i]   -= plank.ENF[i] * ΔT
    @inbounds plank.qNO3[i] -= plank.NR[i] * ΔT
    @inbounds plank.qNH4[i] += plank.NR[i] * ΔT
    @inbounds plank.qNH4[i] += plank.NF[i] * ΔT
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

##### iron allocation
@inline function iron_alloc(par, dpar, qFe, qFePS, qFeNR, qFeNF, p, ΔT)
    # photosynthesis - no-diaz + Crocosphaera
    reg = shape_func_dec(par, p.Imax, 1.0f-2; pow = 4.0f0)
    f_ST2PS_nd = p.k_Fe_ST2PS * reg * qFe * isless(0.0f0, dpar)
    f_PS2ST_nd = p.k_Fe_PS2ST * qFePS * (reg * isless(dpar, 0.0f0) + isequal(0.0f0, par))
    # photosynthesis - trichodesmium
    f_ST2PS_tr = p.k_Fe_ST2PS * qFe * reg * isless(0.0f0, par)
    f_PS2ST_tr = p.k_Fe_PS2ST * qFePS * ((1.0f0 - reg) * isless(0.0f0, par) + isequal(0.0f0, par))
    # photosynthesis
    f_ST2PS = f_ST2PS_nd * (p.is_nr + p.is_croc) + f_ST2PS_tr * p.is_tric
    f_PS2ST = f_PS2ST_nd * (p.is_nr + p.is_croc) + f_PS2ST_tr * p.is_tric

    # nitrate reduction
    f_ST2NR = p.k_Fe_ST2NR * qFe * isequal(0.0f0, par) * p.is_nr 
    f_NR2ST = p.k_Fe_NR2ST * qFeNR * isless(0.0f0, par) * p.is_nr 

    # nitrogen fixation - Crocosphaera watsonii
    f_ST2NF_cr = p.k_Fe_ST2NF * qFe * isequal(0.0f0, par)
    f_NF2ST_cr = p.k_Fe_NF2ST * qFeNF * isless(0.0f0, par)
    # nitrogen fixation - Trichodesmium
    f_ST2NF_tr = p.k_Fe_ST2NF * qFe * (1.0f0 - reg)
    f_NF2ST_tr = p.k_Fe_NF2ST * qFeNF * reg
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

@kernel function calc_iron_fluxes_kernel!(plank, nuts, p, ΔT)
    i = @index(Global)
    @inbounds plank.ST2PS[i], plank.PS2ST[i], plank.ST2NR[i], 
              plank.NR2ST[i], plank.ST2NF[i], plank.NF2ST[i] = 
                    iron_alloc(nuts.par[i], nuts.dpar[i], plank.qFe[i], 
                    plank.qFePS[i], plank.qFeNR[i], plank.qFeNF[i], p, ΔT)
end
function calc_iron_fluxes!(plank, nuts, p, ΔT, arch::Architecture)
    kernel! = calc_iron_fluxes_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p, ΔT)
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
