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
@inline function calc_NH4_uptake(NH4, T, CH, qNH4, Bm, p, ac)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0f-30, NH4+p.KsatNH4) * tempFunc(T, p) * Bm * ac
    return VNH4 * p.is_nr
end

##### calculate nitrate uptake rate (mmolN/individual/second)
@inline function calc_NO3_uptake(NO3, T, CH, qNO3, Bm, p, ac)
    Qn = qNO3/max(1.0f-30, Bm + CH)
    regQN = shape_func_dec(Qn, p.qNO3max, 1.0f-4, pow = 2.0f0)
    VNO3 = p.VNO3max * regQN * NO3/max(1.0f-30, NO3+p.KsatNO3) * tempFunc(T, p) * Bm * ac
    return VNO3 * p.is_nr
end

##### calculate phosphate uptake rate (mmolN/individual/second)
@inline function calc_P_uptake(PO4, T, CH, qP, Bm, p, ac)
    Qp = (qP + Bm * p.R_PC)/max(1.0f-30, Bm + CH)
    regQP = shape_func_dec(Qp, p.qPmax, 1.0f-4, pow = 2.0f0)
    VPO4 = p.VPO4max * regQP * PO4/max(1.0f-30, PO4+p.KsatPO4) * tempFunc(T, p) * Bm * ac
    return VPO4
end

##### calculate iron uptake rate (mmolFe/individual/second)
@inline function calc_Fe_uptake(DFe, T, qFe, qFePS, qFeNR, qFeNF, Bm, CH, Sz, p, ac)
    Qfe = (qFe + qFePS + qFeNF + qFeNR)/max(1.0f-30, Bm + CH)
    regQFe = shape_func_dec(Qfe, p.qFemax, 1.0f-4, pow = 2.0f0)
    SA = 4.0f0 * π * (p.Rad * Sz^(1.0f0/3.0f0))^2.0f0 * 1.0f-12 * p.Nsuper # m²
    VFe = p.KSAFe * SA * DFe * regQFe * tempFunc(T, p) * ac
    return VFe
end

@kernel function calc_NPFe_uptake_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VNH4[i] = calc_NH4_uptake(trs.NH4[i], trs.T[i], plank.CH[i], plank.qNH4[i],
                                              plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VNO3[i] = calc_NO3_uptake(trs.NO3[i], trs.T[i], plank.CH[i], plank.qNO3[i],
                                              plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VPO4[i] = calc_P_uptake(trs.PO4[i], trs.T[i], plank.CH[i], plank.qP[i],
                                            plank.Bm[i], p, plank.ac[i])
    @inbounds plank.VFe[i]  = calc_Fe_uptake(trs.DFe[i], trs.T[i], plank.qFe[i], plank.qFePS[i], 
                                             plank.qFeNR[i], plank.qFeNF[i], plank.Bm[i], plank.CH[i], 
                                             plank.Sz[i], p, plank.ac[i])
end
function calc_NPFe_uptake!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_NPFe_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
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

##### calculate carbon fixation rate (mmolC/individual/second)
@inline function calc_CF(PS, CH, Bm, T, p, ac)
    Qc = CH/max(1.0f-30, Bm + CH)
    regQC = shape_func_dec(Qc, p.CHmax, 1.0f-4, pow = 2.0f0)
    CF = p.k_cf * regQC * tempFunc_CF(T, p) * Bm * ac
    CFt = min(CF, PS / p.e_cf)
    return CFt
end
@kernel function calc_carbon_fixation_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.CF[i] = calc_CF(plank.PS[i], plank.CH[i], 
                                    plank.Bm[i], trs.T[i], p, plank.ac[i])
end
function calc_carbon_fixation!(plank, trs, p, arch::Architecture)
    kernel! = calc_carbon_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate potential nitrate reduction (mmolN/individual/second)
@inline function calc_NR(qNO3, qNH4, qFeNR, Bm, CH, T, p, ac, ΔT)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn , p.qNH4max, 1.0f-4, pow = 2.0f0)
    Qfe_NR = qFeNR / max(1.0f-30, Bm + CH)
    Ksat = Qfe_NR / max(1.0f-30, Qfe_NR + p.KfeNR)
    NR = p.k_nr * reg * qNO3 * Ksat * tempFunc(T, p) * ac
    NR = min(NR, qNO3/ΔT) # double check qNO3 are not over consumed
    return NR * p.is_nr
end
@kernel function calc_NO3_reduction_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.NR[i] = calc_NR(plank.qNO3[i], plank.qNH4[i], plank.qFeNR[i], 
                                    plank.Bm[i], plank.CH[i], trs.T[i], p, plank.ac[i], ΔT)
end
function calc_NO3_reduction!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_NO3_reduction_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate potential nitrogen fixation (mmolN/individual/second)
@inline function calc_NF(qNH4, qFeNF, Bm, CH, T, p, ac)
    Qn = qNH4/max(1.0f-30, Bm + CH)
    reg = shape_func_dec(Qn, p.qNH4max, 1.0f-4, pow = 2.0f0)
    Qfe_NF = qFeNF / max(1.0f-30, Bm + CH)
    Ksat = Qfe_NF / max(1.0f-30, Qfe_NF + p.KfeNF)
    NF = p.k_nf * reg * Ksat * tempFunc(T, p) * Bm * ac
    return NF * (p.is_croc + p.is_tric)
end
@kernel function calc_nitrogen_fixation_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.NF[i] = calc_NF(plank.qNH4[i], plank.qFeNF[i],
                                    plank.Bm[i], plank.CH[i],
                                    trs.T[i], p, plank.ac[i])
end
function calc_nitrogen_fixation!(plank, trs, p, arch::Architecture)
    kernel! = calc_nitrogen_fixation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

@inline function calc_OPS(PS, CF, NR, NF, p)
    # Maximum NADPH required by CF, NF and NR
    NADPH_req = CF * p.re_cf + NF * p.re_nf + NR * p.re_nr
    # Maximum NADPH provided by PS => minimum O₂ production
    OPS = min(PS * p.re_ps, NADPH_req) / p.re_ps * p.o_ps
    return OPS
end
@kernel function calc_O2_production_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.OPS[i] = calc_OPS(plank.PS[i], plank.CF[i],
                                      plank.NR[i], plank.NF[i], p)
end
function calc_O2_production!(plank, p, arch::Architecture)
    kernel! = calc_O2_production_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### passive diffusion of O₂ through cell memberane (positive value means uptake)
@inline function calc_O2_diffusion(O2, qO2, OPS, Sz, p, ΔT, ac)
    Rad = p.Rad * 1.0f-6 # cell radius: m
    Vcell = 4.0f0/3.0f0 * π * (Rad * Sz^(1.0f0/3.0f0))^3.0f0 * p.Nsuper # m³
    O2_incell = (qO2 + OPS * ΔT) / Vcell # concentration mmolO2/m³
    ΔO2 = O2 - O2_incell
    reg = 1.0f0 + shape_func_inc(O2_incell, p.qO2diff, 1.0f-2, pow = 2.0f0)
    VO2 = 4.0f0 * π * Rad * Sz^(1.0f0/3.0f0) * p.k_O2 * ΔO2 * reg * p.Nsuper * ac
    VO2t= VO2 / max(abs(VO2), 1.0f-30) * min(abs(VO2), abs(ΔO2 * Vcell / ΔT))
    return VO2t
end
@kernel function calc_O2_diffusion_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.VO2[i] = calc_O2_diffusion(trs.O2[i], plank.qO2[i], plank.OPS[i],
                                               plank.Sz[i], p, ΔT, plank.ac[i])
end
function calc_O2_diffusion!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_O2_diffusion_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### calculate potential maximum respiration (mmolC/individual/second) 
@inline function calc_respir(PS, CH, qO2, T, p, ac, ΔT)
    RSo = CH * p.k_rs * tempFunc(T, p) * ac
    RSe = CH * p.k_rs * tempFunc(T, p) * ac * isequal(PS, 0.0f0)
    RSo = min(RSo, 0.8f0*CH/ΔT) # double check CH is not over consumed
    RSe = min(RSe, 0.2f0*CH/ΔT) # double check CH is not over consumed
    return RSo, RSe
end
@kernel function calc_respiration_kernel!(plank, trs, p, ΔT)
    i = @index(Global)
    @inbounds plank.RSo[i], plank.RSe[i] = 
                calc_respir(plank.PS[i], plank.CH[i], plank.qO2[i],
                            trs.T[i], p, plank.ac[i], ΔT)
end
function calc_repiration!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_respiration_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### energy allocation
@inline function energy_redox_alloc(PS, OPS, VO2, RSo, RSe, CF, NF, NR, 
                                    ATP, NADPH, qO2, Sz, p, ΔT)
    # Total O₂ available for respiration at current time-step
    O2t = (OPS + VO2) * ΔT + qO2
    RSo_max = min(RSo, O2t/ΔT*0.99f0)
    # ATP and NADPH left after CF from PS
    PSReLeft = max(0.0f0, OPS / p.o_ps * p.re_ps - CF * p.re_cf)
    PSEnLeft = max(0.0f0, PS - CF * p.e_cf)
    # Total ATP and NADPH available
    ReSup_N = PSReLeft + RSe * p.re_rs
    EnSup_N = PSEnLeft + RSo_max * p.e_rs
    # NR and potential NF (without O₂ regulation) based on ATP and NADPH supply
    NRt = min(NR, ReSup_N / p.re_nr)
    NFt = min(NF, ReSup_N / p.re_nf, EnSup_N / p.e_nf)
    # Minimum RSo based on ATP and NDAPH requirement of NF without O₂ regulation
    RSot = min(RSo_max, max(0.0f0, NFt * p.e_nf - PSEnLeft) / p.e_rs)
    # Now, add O₂ regulation on NF, extra energy goes to RSP
    Vcell = 4.0f0/3.0f0 * π * (p.Rad * 1.0f-6 * Sz^(1.0f0/3.0f0))^3.0f0 * p.Nsuper # m³
    O2_incell = (O2t - RSot * ΔT) / Vcell # concentration mmolO2/m³
    regO2 = shape_func_dec(O2_incell, p.qO2nf, 1.0f-2, pow = 2.0f0)
    NFtO2 = NFt * regO2
    # Update RSo and calculate RSP
    RSot2 = min(RSot, max(0.0f0, NFtO2 * p.e_nf - PSEnLeft) / p.e_rs)
    RSP = RSot - RSot2
    # Minimum RSo based on ATP and NDAPH requirement
    RSet = min(RSe, max(0.0f0, NRt * p.re_nr + NFtO2 * p.re_nf - PSReLeft) / p.re_rs)
    # Extra ATP and NADPH
    ATP += ΔT * (PSEnLeft + RSot * p.e_rs - NFt * p.e_nf)
    NADPH += ΔT * (PSReLeft + RSet * p.re_rs - NFt * p.re_nf - NRt * p.re_nr)

    return RSot2, RSet, RSP, NFtO2, NRt, ATP, NADPH
end
@kernel function energy_redox_allocation_kernel!(plank, p, ΔT)
    i = @index(Global)
    plank.RSo[i], plank.RSe[i], plank.RSP[i], plank.NF[i], 
    plank.NR[i], plank.ATP[i], plank.NADPH[i] = 
        energy_redox_alloc(plank.PS[i], plank.OPS[i], plank.VO2[i], 
                           plank.RSo[i], plank.RSe[i],plank.CF[i], 
                           plank.NF[i], plank.NR[i], plank.ATP[i], 
                           plank.NADPH[i], plank.qO2[i], plank.Sz[i], p, ΔT)
end
function energy_redox_allocation!(plank, p, ΔT, arch::Architecture)
    kernel! = energy_redox_allocation_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
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
@inline function iron_alloc(par, dpar, qFe, qFePS, qFeNR, qFeNF, qO2, Sz, tdark, p, ΔT)
    # photosynthesis
    f_ST2PS = p.k_Fe_ST2PS * qFe * isless(0.0f0, dpar)
    f_PS2ST = p.k_Fe_PS2ST * qFePS * (isless(dpar, -10.0f0) + isequal(0.0f0, par))

    # nitrate reduction
    if p.is_nr == 1.0f0
        f_ST2NR = p.k_Fe_ST2NR * qFe
        f_NR2ST = p.k_Fe_NR2ST * qFeNR * (isequal(0.0f0, par) * isless(p.NF_clock, tdark))
    else
        f_ST2NR = 0.0f0
        f_NR2ST = 0.0f0
    end

    # nitrogen fixation
    Rad = p.Rad * 1.0f-6 # cell radius: m
    Vcell = 4.0f0/3.0f0 * π * (Rad * Sz^(1.0f0/3.0f0))^3.0f0 * p.Nsuper # m³
    O2_incell = qO2/Vcell # concentration mmolO2/m³
    regO2 = max(0.0f0, O2_incell / (O2_incell + p.qO2nf))
    if p.is_croc == 1.0f0 # nitrogen fixation - Crocosphaera watsonii
        f_ST2NF = p.k_Fe_ST2NF * qFe * isequal(0.0f0, par) * isless(tdark, p.NF_clock)
        f_NF2ST = p.k_Fe_NF2ST * qFeNF * (isequal(0.0f0, par) * 
                    isless(p.NF_clock, tdark) + isless(0.0f0, par) + regO2)
    elseif p.is_tric ==1.0f0 # nitrogen fixation - Trichodesmium
        f_ST2NF = p.k_Fe_ST2NF * qFe * isless(0.0f0, dpar)
        f_NF2ST = p.k_Fe_NF2ST * qFeNF * 
                    (isless(dpar, 0.0f0) + isequal(0.0f0, par) + regO2)
    else
        f_ST2NF = 0.0f0
        f_NF2ST = 0.0f0
    end

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
                               plank.qFePS[i], plank.qFeNR[i], 
                               plank.qFeNF[i], plank.qO2[i], plank.Sz[i],
                               plank.tdark[i], p, ΔT)
end
function calc_iron_fluxes!(plank, trs, p, ΔT, arch::Architecture)
    kernel! = calc_iron_fluxes_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p, ΔT)
    return nothing
end

##### update C, N, P, Fe quotas, biomass, Chla
@kernel function update_states_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.qP[i]   += ΔT * (plank.VPO4[i] - plank.BS[i] * p.R_PC)
    @inbounds plank.qNO3[i] += ΔT * (plank.VNO3[i] - plank.NR[i])
    @inbounds plank.qNH4[i] += ΔT * (plank.VNH4[i] + plank.NR[i] +
                                     plank.NF[i]   - plank.BS[i] * p.R_NC)
    @inbounds plank.qO2[i]  += ΔT * (plank.OPS[i]  + plank.VO2[i] -
                                     plank.RSo[i]  - plank.RSP[i])
    @inbounds plank.CH[i]   += ΔT * (plank.CF[i]   - plank.RSo[i] - 
                                     plank.RSe[i]  - plank.RSP[i] -
                                     plank.BS[i])
    @inbounds plank.Bm[i]   += ΔT *  plank.BS[i]
    @inbounds plank.qFe[i]  += ΔT * (plank.PS2ST[i] - plank.ST2PS[i] +
                                     plank.NR2ST[i] - plank.ST2NR[i] +
                                     plank.NF2ST[i] - plank.ST2NF[i] +
                                     plank.VFe[i])
    @inbounds plank.qFePS[i]+= ΔT * (plank.ST2PS[i] - plank.PS2ST[i])
    @inbounds plank.qFeNR[i]+= ΔT * (plank.ST2NR[i] - plank.NR2ST[i])
    @inbounds plank.qFeNF[i]+= ΔT * (plank.ST2NF[i] - plank.NF2ST[i])

    @inbounds plank.Chl[i]  += ΔT *  plank.BS[i] * p.R_NC * plank.ρChl[i]
    @inbounds plank.age[i]  += ΔT /  3600.0f0 * plank.ac[i]
end
function update_states!(plank, p, ΔT, arch::Architecture)
    kernel! = update_states_kernel!(device(arch), 256, (size(plank.ac,1)))
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