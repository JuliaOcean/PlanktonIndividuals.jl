##### set up the operating array(cuarray) for plankton physiology processes
#= index
1  2  3  4       5     6   7   8   9   10   11  12   13   14   15   16   17
x  y  z  size_i  size  Bm  Cq  Nq  Pq  chl  sp  gen  age  NH4  NO3  PO4  DOC

18  19        20  21    22    23    24    25    26    27  28   29   30    31
αI  TempFunc  PS  VDOC  VNH4  VNO3  VPO4  ρchl  resp  BS  exu  grz  mort  dvid
=#
function phyt_op_array_setup(phytos, arch::Architecture)
    total_num = size(phytos, 1)
    op_array = zeros(total_num, 31) |> array_type(arch)

    op_array[:, 1:13] .= phytos[:, 1:13]
    return op_array
end

##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(op_array, NH4, NO3, PO4, DOC, par, temp, g::Grids, α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    sp = op_array[i,11]
    op_array[i,14] = max(0.0, NH4[xi, yi, zi])
    op_array[i,15] = max(0.0, NO3[xi, yi, zi])
    op_array[i,16] = max(0.0, PO4[xi, yi, zi])
    op_array[i,17] = max(0.0, DOC[xi, yi, zi])
    op_array[i,18] = α[sp] * par[xi, yi, zi] * Φ[sp]
    op_array[i,19] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15) - 1.0 / Tempref))) * TempCoeff
end
function find_NPT!(op_array, arch::Architecture, NH4, NO3, PO4, DOC, par, temp, g::Grids,
                   α, Φ, TempAe, Tempref, TempCoeff)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(op_array,1),))

    event = kernel!(op_array, NH4, NO3, PO4, DOC, par, temp, g,
                    α, Φ, TempAe, Tempref, TempCoeff)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@kernel function calc_PS_kernel!(op_array, PCmax, PC_b)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    PCm = PCmax[sp] * op_array[i,19] * op_array[i,5]^PC_b[sp]
    op_array[i,20] = PCm * (1.0 - exp(-op_array[i,18] * op_array[i,10] / op_array[i,6] / PCm)) * op_array[i,6]
end
function calc_PS!(op_array, arch::Architecture, PCmax, PC_b)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, PCmax, PC_b)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
@kernel function calc_VDOC_kernel!(op_array, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    i = @index(Global, Linear)
    sp = op_array[i,11]

    Qc = op_array[i,7] / (op_array[i,6] + op_array[i,7])
    reg= max(0.0, min(1.0, (Cqmax[sp] - Qc) / (Cqmax[sp] - Cqmin[sp])))

    ksat = op_array[i,17] / (op_array[i,17] + KsatDOC[sp])

    op_array[i,21] = VDOCmax[sp] * op_array[i,5]^VDOC_b[sp] * ksat * reg * op_array[i,19] * op_array[i,6]
    op_array[i,21] = min(op_array[i,17] * g.V / 10.0 / ΔT, op_array[i,21])
end
function calc_VDOC!(op_array, arch::Architecture, DOC, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    kernel! = calc_VDOC_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, g, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    wait(device(arch), event)
    return nothing
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
@kernel function calc_VN_kernel!(op_array, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    i = @index(Global, Linear)
    sp = op_array[i,11]

    Qn = (op_array[i,8] + op_array[i,6] * R_NC) / (op_array[i,6] + op_array[i,7])
    reg= max(0.0, min(1.0, (Nqmax[sp] - Qn) / (Nqmax[sp] - Nqmin[sp])))

    ksatNH = op_array[i,14] / (op_array[i,14] + KsatNH4[sp])
    ksatNO = op_array[i,15] / (op_array[i,15] + KsatNO3[sp])

    op_array[i,22] = VNH4max[sp] * op_array[i,5]^VN_b[sp] * ksatNH * reg * op_array[i,19] * op_array[i,6]
    op_array[i,22] = min(op_array[i,14] * g.V / 10.0 / ΔT, op_array[i,22])

    op_array[i,23] = VNO3max[sp] * op_array[i,5]^VN_b[sp] * ksatNO * reg * op_array[i,19] * op_array[i,6]
    op_array[i,23] = min(op_array[i,15] * g.V / 10.0 / ΔT, op_array[i,23])
end
function calc_VN!(op_array, arch::Architecture, g::Grids, ΔT,
                  Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    kernel! = calc_VN_kernel!(device(arch), 256, (size(op_array,1),))

    event = kernel!(op_array, g, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    wait(device(arch), event)
    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
@kernel function calc_VP_kernel!(op_array, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    i = @index(Global, Linear)
    sp = op_array[i,11]

    Qp = (op_array[i,9] + op_array[i,6] * R_PC) / (op_array[i,6] + op_array[i,7])
    reg= max(0.0, min(1.0, (Pqmax[sp] - Qp) / (Pqmax[sp] - Pqmin[sp])))

    ksat = op_array[i,16] / (op_array[i,16] + KsatPO4[sp])

    op_array[i,24] = VPO4max[sp] * op_array[i,5]^VP_b[sp] * ksat * reg * op_array[i,19] * op_array[i,6]
    op_array[i,24] = min(op_array[i,16] * g.V / 10.0 / ΔT, op_array[i,24])
end
function calc_VP!(op_array, arch::Architecture, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    kernel! = calc_VP_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, g, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(op_array, Chl2N)
    i = @index(Global, Linear)
    op_array[i,25] = op_array[i,18] < 0 ? 0.0 :
        op_array[i,20] / op_array[i,6] * Chl2N / (op_array[i,18] * op_array[i,10] / op_array[i,6])
end
function calc_ρchl!(op_array, arch::Architecture, Chl2N)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(op_array, respir_a, respir_b)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    op_array[i,26] = respir_a[sp] * op_array[i,5]^respir_b[sp] * op_array[i,6] * op_array[i,19]
end
function calc_respir!(op_array, arch::Architecture, respir_a, respir_b)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, respir_a, respir_b)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas
@kernel function update_quotas_kernel!(op_array, R_NC, R_PC, ΔT)
    i = @index(Global, Linear)
    op_array[i,7] = op_array[i,7] + ΔT * (op_array[i,20] + op_array[i,21] - op_array[i,26])
    if op_array[i,7] ≤ 0.0  # if C reserve is not enough for respiration
        exceed = 0.0 - op_array[i,7]
        op_array[i,7] = 0.0
        op_array[i,8] = op_array[i,8] + exceed * R_NC # return N from function pool to N reserve
        op_array[i,9] = op_array[i,9] + exceed * R_PC # return P from function pool to P reserve
    end
    op_array[i,8] = op_array[i,8] + ΔT * (op_array[i,22] + op_array[i,23])
    op_array[i,9] = op_array[i,9] + ΔT *  op_array[i,24]
end
function update_quotas!(op_array, arch::Architecture, R_NC, R_PC, ΔT)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, R_NC, R_PC, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(op_array, k_mtb, b_k_mtb, R_NC, R_PC)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    BS_Cmax = k_mtb[sp] * op_array[i,5]^b_k_mtb[sp] * op_array[i,7]
    op_array[i,27] = min(op_array[i,7], op_array[i,8]/R_NC, op_array[i,9]/R_PC) * k_mtb[sp] * op_array[i,5]^b_k_mtb[sp]
    op_array[i,28] = BS_Cmax - op_array[i,27]
end
function calc_BS!(op_array, arch::Architecture, k_mtb, b_k_mtb, R_NC, R_PC)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, k_mtb, b_k_mtb, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
@kernel function update_biomass_kernel!(op_array, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    op_array[i,6] = op_array[i,6]  + ΔT *  op_array[i,27]
    op_array[i,7] = op_array[i,7]  - ΔT * (op_array[i,27] + op_array[i,28])
    op_array[i,8] = op_array[i,8]  - ΔT *  op_array[i,27] * R_NC
    op_array[i,9] = op_array[i,9]  - ΔT *  op_array[i,27] * R_PC
    op_array[i,10]= op_array[i,10] + ΔT *  op_array[i,27] * R_NC * op_array[i,25]
    op_array[i,13]= op_array[i,13] + ΔT /3600
    op_array[i,5] =(op_array[i,6]  + op_array[i,7]) / P_Cquota[sp] / P_Nsuper

end
function update_biomass!(op_array, arch::Architecture, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(op_array,1),))
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
        reg_graz = 1.0 / Grz_P * max(0.15, 1 - abs(op_array[i,3]) / Grz_stp)
    end
    op_array[i,29] = rand(Bernoulli(reg_graz))
end
function calc_graz!(op_array, arch::Architecture, Grz_P, Grz_stp)
    kernel! = calc_graz_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, Grz_P, Grz_stp)
    wait(device(arch), event)
    return nothing
end


##### calculate probability of natural death
@kernel function calc_mort_kernel!(op_array, mort_reg, mort_P)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    reg = mort_P[sp] * (tanh(6.0 * (mort_reg[sp] - op_array[i,5])) + 1)
    op_array[30,i] = rand(Bernoulli(reg))
end
function calc_mort!(op_array, arch::Architecture, mort_reg, mort_P)
    kernel! = calc_mort_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, mort_reg, mort_P)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of sizer-like cell division
@kernel function calc_dvid_kernel!(op_array, dvid_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    i = @index(Global, Linear)
    sp = op_array[i,11]
    if op_array[i,6] ≥ 2 * P_Cquota[sp] * P_Nsuper
        if dvid_type[sp] == 1
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[i,5] - dvid_reg[sp]))+1)
            op_array[i,31] = rand(Bernoulli(reg))
        elseif dvid_type == 2
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[i,5] - op_array[i,4] - dvid_reg[sp]))+1)
            op_array[i,31] = rand(Bernoulli(reg))
        elseif dvid_type == 3
            reg = dvid_P[sp] * (tanh(dvid_stp[sp] * (op_array[i,13] - dvid_reg[sp]))+1)
            op_array[i,31] = rand(Bernoulli(reg))
        elseif dvid_type == 4
            cirT = t % 86400 ÷ 3600
            reg  = dvid_P[sp] * (tanh(dvid_stp[sp] * (cirT - dvid_reg[sp]))+1)
            op_array[i,31] = rand(Bernoulli(reg))
        elseif dvid_type == 5
            cirT = t % 86400 ÷ 3600
            regT = dvid_stp[sp][1] * (cirT - dvid_reg[sp][1])
            regS = dvid_stp[sp][2] * (op_array[i,5] - dvid_reg[sp][2])
            reg  = dvid_P[sp] * (tanh(regS) + 1) * (tanh(regT) + 1)
            op_array[i,31] = rand(Bernoulli(reg))
        else
            throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
        end
    end
end
function calc_dvid!(op_array, arch::Architecture, dvid_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    kernel! = calc_dvid_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, dvie_type, dvid_stp, dvid_P, dvid_reg, P_Cquota, P_Nsuper)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(op_array, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (op_array[i,6] + op_array[i,7]) * lossFracC
    POC_con[xi, yi, zi] = POC_con[xi, yi, zi] + (op_array[i,6] + op_array[i,7]) * (1.0 - lossFracC)
    DON_con[xi, yi, zi] = DON_con[xi, yi, zi] + (op_array[i,6] * R_NC + op_array[i,8]) * lossFracN
    PON_con[xi, yi, zi] = PON_con[xi, yi, zi] + (op_array[i,6] * R_NC + op_array[i,8]) * (1.0 - lossFracN)
    DOP_con[xi, yi, zi] = DOP_con[xi, yi, zi] + (op_array[i,6] * R_PC + op_array[i,9]) * lossFracP
    POP_con[xi, yi, zi] = POP_con[xi, yi, zi] + (op_array[i,6] * R_PC + op_array[i,9]) * (1.0 - lossFracP)
end
function calc_loss!(op_array, arch::Architecture, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_P)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(op_array, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    xi = find_xF_ind(op_array[i,1], g) |> Int
    yi = find_yF_ind(op_array[i,2], g) |> Int
    zi = find_zF_ind(op_array[i,3], g) |> Int
    DIC_con[xi, yi, zi] = DIC_con[xi, yi, zi] + (op_array[i,26] - op_array[20,i]) * ΔT
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (op_array[i,30] - op_array[21,i]) * ΔT
    NH4_con[xi, yi, zi] = NH4_con[xi, yi, zi] -  op_array[i,22] * ΔT
    NO3_con[xi, yi, zi] = NO3_con[xi, yi, zi] -  op_array[i,23] * ΔT
    PO4_con[xi, yi, zi] = PO4_con[xi, yi, zi] -  op_array[i,24] * ΔT
end
function calc_consume!(op_array, arch::Architecture, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(op_array,1),))
    event = kernel!(op_array, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g, ΔT)
    wait(device(arch), event)
    return nothing
end
