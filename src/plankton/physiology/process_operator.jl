##### calculate αI for each individual
@kernel function calc_αI_kernel!(op_array, par, g::Grids, α, Φ)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    sp = op_array[11,i]
    op_array[14,i] = α[sp] * par[xi, yi, zi] * Φ[sp]
end
function calc_αI!(op_array, arch::Architecture, par, g::Grids, α, Φ)
    kernel! = calc_αI_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, par, α, Φ, g)
    wait(device(arch), event)
    return nothing
end

##### calculate the limitation of temperature
@kernel function calc_tempfunc_kernel!(op_array, temp, g::Grids, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    op_array[15,i] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15) - 1.0 / Tempref))) * TempCoeff
end
function calc_tempfunc!(op_array, arch::Architecture, temp, g::Grids, TempAe, Tempref, TempCoeff)
    kernel! = calc_tempfunc_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, temp, TempAe, Tempref, TempCoeff, g)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@kernel function calc_PS_kernel!(op_array, PCmax, PC_b)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    PCm = PCmax[sp] * op_array[15,i] * op_array[5,i]^PC_b[sp]
    op_array[16,i] = PCm * (1.0 - exp(-op_array[14,i] * op_array[10,i] / op_array[6,i] / PCm)) * op_array[6,i]
end
function calc_PS!(op_array, arch::Architecture, PCmax, PC_b)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, PCmax, PC_b)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
@kernel function calc_VDOC_kernel!(op_array, DOC, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    sp = op_array[11,i]
    Qc = op_array[7,i] / (op_array[6,i] + op_array[7,i])
    reg= max(0.0, min(1.0, (Cqmax[sp] - Qc) / (Cqmax[sp] - Cqmin[sp])))
    ksat = DOC[xi, yi, zi] / (DOC[xi, yi, zi] + KsatDOC[sp])
    op_array[17,i] = VDOCmax[sp] * op_array[5,i]^VDOC_b[sp] * ksat * reg * op_array[15,i] * op_array[6,i]
    op_array[17,i] = min(DOC[xi, yi, zi] * g.V / 10.0 / ΔT, op_array[17,i])
end
function calc_VDOC!(op_array, arch::Architecture, DOC, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    kernel! = calc_VDOC_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, DOC, g, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    wait(device(arch), event)
    return nothing
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
@kernel function calc_VN_kernel!(op_array, NH4, NO3, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    sp = op_array[11,i]
    Qn = (op_array[8,i] + op_array[6,i] * R_NC) / (op_array[6,i] + op_array[7,i])
    reg= max(0.0, min(1.0, (Nqmax[sp] - Qn) / (Nqmax[sp] - Nqmin[sp])))
    ksatNH = NH4[xi, yi, zi] / (NH4[xi, yi, zi] + KsatNH4[sp])
    ksatNO = NO3[xi, yi, zi] / (NO3[xi, yi, zi] + KsatNO3[sp])
    op_array[18,i] = VNH4max[sp] * op_array[5,i]^VN_b[sp] * ksatNH * reg * op_array[15,i] * op_array[6,i]
    op_array[18,i] = min(NH4[xi, yi, zi] * g.V / 10.0 / ΔT, op_array[18,i])
    op_array[19,i] = VNO3max[sp] * op_array[5,i]^VN_b[sp] * ksatNO * reg * op_array[15,i] * op_array[6,i]
    op_array[19,i] = min(NO3[xi, yi, zi] * g.V / 10.0 / ΔT, op_array[19,i])
end
function calc_VN!(op_array, arch::Architecture, NH4, NO3, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    kernel! = calc_VN_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, NH4, NO3, g, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    wait(device(arch), event)
    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
@kernel function calc_VP_kernel!(op_array, PO4, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    sp = op_array[11,i]
    Qp = (op_array[9,i] + op_array[6,i] * R_PC) / (op_array[6,i] + op_array[7,i])
    reg= max(0.0, min(1.0, (Pqmax[sp] - Qp) / (Pqmax[sp] - Pqmin[sp])))
    ksat = PO4[xi, yi, zi] / (PO4[xi, yi, zi] + KsatPO4[sp])
    op_array[20,i] = VPO4max[sp] * op_array[5,i]^VP_b[sp] * ksat * reg * op_array[15,i] * op_array[6,i]
    op_array[20,i] = min(PO4[xi, yi, zi] * g.V / 10.0 / ΔT, op_array[20,i])
end
function calc_VP!(op_array, arch::Architecture, PO4, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    kernel! = calc_VP_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, PO4, g, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(op_array, Chl2N)
    i = @index(Global, Linear)
    if op_array[14,i] > 0
        op_array[21,i] = op_array[16,i] / op_array[6,i] * Chl2N / (op_array[14,i] * op_array[10,i] / op_array[6,i])
    else
        op_array[21,i] = 0.0
    end
end
function calc_ρchl!(op_array, arch::Architecture, Chl2N)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(op_array, respir_a, respir_b)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    op_array[22,i] = respir_a[sp] * op_array[5,i]^respir_b[sp] * op_array[6,i] * op_array[15,i]
end
function calc_respir!(op_array, arch::Architecture, respir_a, respir_b)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, respir_a, respir_b)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas
@kernel function update_quotas_kernel!(op_array, R_NC, R_PC, ΔT)
    i = @index(Global, Linear)
    op_array[7,i] = op_array[7,i] + ΔT * (op_array[16,i] + op_array[17,i] - op_array[22,i])
    if op_array[7,i] ≤ 0.0  # if C reserve is not enough for respiration
        exceed = 0.0 - op_array[7,i]
        op_array[7,i] = 0.0
        op_array[8,i] = op_array[8,i] + exceed * R_NC # return N from function pool to N reserve
        op_array[9,i] = op_array[9,i] + exceed * R_PC # return P from function pool to P reserve
    end
    op_array[8,i] = op_array[8,i] + ΔT * (op_array[18,i] + op_array[19,i])
    op_array[9,i] = op_array[9,i] + ΔT *  op_array[20,i]
end
function update_quotas!(op_array, arch::Architecture, R_NC, R_PC, ΔT)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, R_NC, R_PC, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(op_array, k_mtb, b_k_mtb, R_NC, R_PC)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    BS_Cmax = k_mtb[sp] * op_array[5,i]^b_k_mtb[sp] * op_array[7,i]
    op_array[23,i] = min(op_array[7,i], op_array[8,i]/R_NC, op_array[9,i]/R_PC) * k_mtb[sp] * op_array[5,i]^b_k_mtb[sp]
    op_array[24,i] = BS_Cmax - op_array[23,i]
end
function calc_BS!(op_array, arch::Architecture, k_mtb, b_k_mtb, R_NC, R_PC)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, k_mtb, b_k_mtb, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
@kernel function update_biomass_kernel!(op_array, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    op_array[6,i] = op_array[6,i]  + ΔT *  op_array[23,i]
    op_array[7,i] = op_array[7,i]  - ΔT * (op_array[23,i] + op_array[24,i])
    op_array[8,i] = op_array[8,i]  - ΔT *  op_array[23,i] * R_NC
    op_array[9,i] = op_array[9,i]  - ΔT *  op_array[23,i] * R_PC
    op_array[10,i]= op_array[10,i] + ΔT *  op_array[23,i] * R_NC * op_array[21,i]
    op_array[13,i]= op_array[13,i] + ΔT /3600
    op_array[5,i] =(op_array[6,i]  + op_array[7,i]) / P_Cquota[sp] / P_Nsuper

end
function update_biomass!(op_array, arch::Architecture, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of grazing
@kernel function calc_graz_kernel!(op_array, Grz_P, Grz_stp)
    i = @index(Global, Linear)
    if Grz_stp == 0
        reg_graz = 1.0 / Grz_P
    else
        reg_graz = 1.0 / Grz_P * max(0.15, 1 - abs(op_array[3,i]) / Grz_stp)
    end
    op_array[25,i] = rand(Bernoulli(reg_graz))
end
function calc_graz!(op_array, arch::Architecture, Grz_P, Grz_stp)
    kernel! = calc_graz_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, Grz_P, Grz_stp)
    wait(device(arch), event)
    return nothing
end


##### calculate probability of natural death
@kernel function calc_mort_kernel!(op_array, mort_reg, mort_P)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    reg = mort_P[sp] * (tanh(6.0 * (mort_reg[sp] - op_array[5,i])) + 1)
    op_array[26,i] = rand(Bernoulli(reg))
end
function calc_mort!(op_array, arch::Architecture, mort_reg, mort_P)
    kernel! = calc_mort_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, mort_reg, mort_P)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of sizer-like cell division
@kernel function calc_dvid_kernel!(op_array, dvid_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    i = @index(Global, Linear)
    sp = op_array[11,i]
    if op_array[6,i] ≥ 2 * P_Cquota[sp] * P_Nsuper
        if dvid_type[sp] == 1
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[5,i] - dvid_reg[sp]))+1)
            op_array[27,i] = rand(Bernoulli(reg))
        elseif dvid_type == 2
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[5,i] - op_array[4,i] - dvid_reg[sp]))+1)
            op_array[27,i] = rand(Bernoulli(reg))
        elseif dvid_type == 3
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[13,i] - dvid_reg[sp]))+1)
            op_array[27,i] = rand(Bernoulli(reg))
        elseif dvid_type == 4
            cirT = t % 86400 ÷ 3600
            reg  = dvid_P[sp] * (tanh(dvid_stp[sp] * (cirT - dvid_reg[sp]))+1)
            op_array[27,i] = rand(Bernoulli(reg))
        elseif dvid_type == 5
            cirT = t % 86400 ÷ 3600
            regT = dvid_stp[sp][1] * (cirT - dvid_reg[sp][1])
            regS = dvid_stp[sp][2] * (op_array[5,i] - dvid_reg[sp][2])
            reg  = dvid_P[sp] * (tanh(regS) + 1) * (tanh(regT) + 1)
            op_array[27,i] = rand(Bernoulli(reg))
        else
            throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
        end
    end
end
function calc_dvid!(op_array, arch::Architecture, dvid_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    kernel! = calc_dvid_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, dvie_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(op_array, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (op_array[6,i] + op_array[7,i]) * lossFracC
    POC_con[xi, yi, zi] = POC_con[xi, yi, zi] + (op_array[6,i] + op_array[7,i]) * (1.0 - lossFracC)
    DON_con[xi, yi, zi] = DON_con[xi, yi, zi] + (op_array[6,i] * R_NC + op_array[8,i]) * lossFracN
    PON_con[xi, yi, zi] = PON_con[xi, yi, zi] + (op_array[6,i] * R_NC + op_array[8,i]) * (1.0 - lossFracN)
    DOP_con[xi, yi, zi] = DOP_con[xi, yi, zi] + (op_array[6,i] * R_PC + op_array[9,i]) * lossFracP
    POP_con[xi, yi, zi] = POP_con[xi, yi, zi] + (op_array[6,i] * R_PC + op_array[9,i]) * (1.0 - lossFracP)
end
function calc_loss!(op_array, arch::Architecture, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_P)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(op_array, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[1,i], g)
    yi = find_yF_ind(op_array[2,i], g)
    zi = find_zF_ind(op_array[3,i], g)
    DIC_con[xi, yi, zi] = DIC_con[xi, yi, zi] + (op_array[22,i] - op_array[16,i]) * ΔT
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (op_array[24,i] - op_array[17,i]) * ΔT
    NH4_con[xi, yi, zi] = NH4_con[xi, yi, zi] -  op_array[18,i] * ΔT
    NO3_con[xi, yi, zi] = NO3_con[xi, yi, zi] -  op_array[19,i] * ΔT
    PO4_con[xi, yi, zi] = PO4_con[xi, yi, zi] -  op_array[20,i] * ΔT
end
function calc_consume!(op_array, arch::Architecture, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(op_array,2),))
    event = kernel!(op_array, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g, ΔT)
    wait(device(arch), event)
    return nothing
end
