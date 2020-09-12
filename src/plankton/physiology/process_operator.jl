##### set up the operating array(cuarray) for plankton physiology processes
#= index
1  2  3  4       5     6   7   8   9   10   11  12   13   14  15  16  17   18   19   20
x  y  z  size_i  size  Bm  Cq  Nq  Pq  chl  sp  gen  age  xi  yi  zi  NH4  NO3  PO4  DOC

21  22        23  24    25    26    27    28    29    30  31   32   33    34
αI  TempFunc  PS  VDOC  VNH4  VNO3  VPO4  ρchl  resp  BS  exu  grz  mort  dvid
=#
@kernel function update_phy_ope_kernel!(phytos, ope)
    i = @index(Global, Linear)
    ope[i,1]  = phytos[i,1]
    ope[i,2]  = phytos[i,2]
    ope[i,3]  = phytos[i,3]
    ope[i,4]  = phytos[i,4]
    ope[i,5]  = phytos[i,5]
    ope[i,6]  = phytos[i,6]
    ope[i,7]  = phytos[i,7]
    ope[i,8]  = phytos[i,8]
    ope[i,9]  = phytos[i,9]
    ope[i,10] = phytos[i,10]
    ope[i,11] = phytos[i,11]
    ope[i,12] = phytos[i,12]
    ope[i,13] = phytos[i,13]
end
function update_phy_ope!(phytos, ope, arch::Architecture)
    kernel! = update_phy_ope_kernel!(device(arch), 256, (size(phytos,1),))
    event = kernel!(phytos,ope)
    wait(device(arch), event)
    return nothing
end

##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(ope, inds::AbstractArray{Int64,2}, sp::AbstractArray{Int64,1},
                                  NH4, NO3, PO4, DOC, par, temp, g::Grids,
                                  α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    sp = sp[i]
    ope[i,17] = max(0.0, NH4[xi, yi, zi])
    ope[i,18] = max(0.0, NO3[xi, yi, zi])
    ope[i,19] = max(0.0, PO4[xi, yi, zi])
    ope[i,20] = max(0.0, DOC[xi, yi, zi])
    ope[i,21] = α[sp] * par[xi, yi, zi] * Φ[sp]
    ope[i,22] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15) - 1.0 / Tempref))) * TempCoeff
end
function find_NPT!(ope, inds::AbstractArray{Int64,2}, sp::AbstractArray{Int64,1}, arch::Architecture,
                   NH4, NO3, PO4, DOC, par, temp, g::Grids,
                   α, Φ, TempAe, Tempref, TempCoeff)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(ope,1),))

    event = kernel!(ope, inds, sp,  NH4, NO3, PO4, DOC, par, temp, g,
                    α, Φ, TempAe, Tempref, TempCoeff)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@kernel function calc_PS_kernel!(ope, sp::AbstractArray{Int64,1}, PCmax, PC_b)
    i = @index(Global, Linear)
    sp = sp[i]
    PCm = PCmax[sp] * ope[i,5]^PC_b[sp] * ope[i,22]
    ope[i,23] = PCm * (1.0 - exp(-ope[i,21] * ope[i,10] / ope[i,6] / PCm)) * ope[i,6]
end
function calc_PS!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, PCmax, PC_b)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, PCmax, PC_b)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
@kernel function calc_VDOC_kernel!(ope, sp::AbstractArray{Int64,1}, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    i = @index(Global, Linear)
    sp = sp[i]

    Qc = ope[i,7] / (ope[i,6] + ope[i,7])
    reg= max(0.0, min(1.0, (Cqmax[sp] - Qc) / (Cqmax[sp] - Cqmin[sp])))

    ksat = ope[i,20] / (ope[i,20] + KsatDOC[sp])

    ope[i,24] = VDOCmax[sp] * ope[i,5]^VDOC_b[sp] * ksat * reg * ope[i,22] * ope[i,6]
    ope[i,24] = min(ope[i,20] * g.V / 10.0 / ΔT, ope[i,24])
end
function calc_VDOC!(ope, sp::AbstractArray{Int64,1}, arch::Architecture,
                    g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    kernel! = calc_VDOC_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, g, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    wait(device(arch), event)
    return nothing
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
@kernel function calc_VN_kernel!(ope, sp::AbstractArray{Int64,1}, g::Grids,
                                 ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    i = @index(Global, Linear)
    sp = sp[i]

    Qn = (ope[i,8] + ope[i,6] * R_NC) / (ope[i,6] + ope[i,7])
    reg= max(0.0, min(1.0, (Nqmax[sp] - Qn) / (Nqmax[sp] - Nqmin[sp])))

    ksatNH = ope[i,17] / (ope[i,17] + KsatNH4[sp])
    ksatNO = ope[i,18] / (ope[i,18] + KsatNO3[sp])

    ope[i,25] = VNH4max[sp] * ope[i,5]^VN_b[sp] * ksatNH * reg * ope[i,22] * ope[i,6]
    ope[i,25] = min(ope[i,17] * g.V / 10.0 / ΔT, ope[i,25])

    ope[i,26] = VNO3max[sp] * ope[i,5]^VN_b[sp] * ksatNO * reg * ope[i,22] * ope[i,6]
    ope[i,26] = min(ope[i,18] * g.V / 10.0 / ΔT, ope[i,26])
end
function calc_VN!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, g::Grids, ΔT,
                  Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    kernel! = calc_VN_kernel!(device(arch), 256, (size(ope,1),))

    event = kernel!(ope, sp, g, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    wait(device(arch), event)
    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
@kernel function calc_VP_kernel!(ope, sp::AbstractArray{Int64,1}, g::Grids,
                                 ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    i = @index(Global, Linear)
    sp = sp[i]

    Qp = (ope[i,9] + ope[i,6] * R_PC) / (ope[i,6] + ope[i,7])
    reg= max(0.0, min(1.0, (Pqmax[sp] - Qp) / (Pqmax[sp] - Pqmin[sp])))

    ksat = ope[i,19] / (ope[i,19] + KsatPO4[sp])

    ope[i,27] = VPO4max[sp] * ope[i,5]^VP_b[sp] * ksat * reg * ope[i,22] * ope[i,6]
    ope[i,27] = min(ope[i,19] * g.V / 10.0 / ΔT, ope[i,27])
end
function calc_VP!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, g::Grids,
                  ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    kernel! = calc_VP_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, g, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(ope, Chl2N)
    i = @index(Global, Linear)
    ope[i,28] = ope[i,21] < 0 ? 0.0 :
        ope[i,23] / ope[i,6] * Chl2N / (ope[i,21] * ope[i,10] / ope[i,6])
end
function calc_ρchl!(ope, arch::Architecture, Chl2N)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(ope, sp::AbstractArray{Int64,1}, respir_a, respir_b)
    i = @index(Global, Linear)
    sp = sp[i]
    ope[i,29] = respir_a[sp] * ope[i,5]^respir_b[sp] * ope[i,6] * ope[i,22]
end
function calc_respir!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, respir_a, respir_b)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, respir_a, respir_b)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas
@kernel function update_quotas_kernel!(ope, R_NC, R_PC, ΔT)
    i = @index(Global, Linear)
    ope[i,7] = ope[i,7] + ΔT * (ope[i,23] + ope[i,24] - ope[i,29])
    if ope[i,7] ≤ 0.0  # if C reserve is not enough for respiration
        exceed = 0.0 - ope[i,7]
        ope[i,7] = 0.0
        ope[i,8] = ope[i,8] + exceed * R_NC # return N from function pool to N reserve
        ope[i,9] = ope[i,9] + exceed * R_PC # return P from function pool to P reserve
    end
    ope[i,8] = ope[i,8] + ΔT * (ope[i,25] + ope[i,26])
    ope[i,9] = ope[i,9] + ΔT *  ope[i,27]
end
function update_quotas!(ope, arch::Architecture, R_NC, R_PC, ΔT)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, R_NC, R_PC, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(ope, sp::AbstractArray{Int64,1}, k_mtb, b_k_mtb, R_NC, R_PC)
    i = @index(Global, Linear)
    sp = sp[i]
    BS_Cmax = k_mtb[sp] * ope[i,5]^b_k_mtb[sp] * ope[i,7]
    ope[i,30] = min(ope[i,7], ope[i,8]/R_NC, ope[i,9]/R_PC) * k_mtb[sp] * ope[i,5]^b_k_mtb[sp]
    ope[i,31] = max(0.0, BS_Cmax - ope[i,30])
end
function calc_BS!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, k_mtb, b_k_mtb, R_NC, R_PC)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, k_mtb, b_k_mtb, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
@kernel function update_biomass_kernel!(ope, sp::AbstractArray{Int64,1}, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    i = @index(Global, Linear)
    sp = sp[i]
    ope[i,6] = ope[i,6]  + ΔT *  ope[i,30]
    ope[i,7] = ope[i,7]  - ΔT * (ope[i,30] + ope[i,31])
    ope[i,8] = ope[i,8]  - ΔT *  ope[i,30] * R_NC
    ope[i,9] = ope[i,9]  - ΔT *  ope[i,30] * R_PC
    ope[i,10]= ope[i,10] + ΔT *  ope[i,30] * R_NC * ope[i,28]
    ope[i,13]= ope[i,13] + ΔT /3600
    ope[i,5] =(ope[i,6]  + ope[i,7]) / P_Cquota[sp] / P_Nsuper

end
function update_biomass!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of grazing, mortality and cell division
@kernel function probabilities_kernel!(ope, sp::AbstractArray{Int64,1}, Grz_P, Grz_stp,
                                       mort_reg, mort_P, dvid_type, dvid_stp, dvid_stp2,
                                       dvid_P, dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    i = @index(Global, Linear)
    sp = sp[i]

    if t%300 == 1
        ##### grazing
        if Grz_P == 0
            ope[i,32] = 0.0
        else
            reg_grz = Grz_stp == 0 ? 1.0 / Grz_P : 1.0 / Grz_P * max(0.15, 1 - abs(ope[i,3]) / Grz_stp)
            # reg = 1.0 / Grz_P
            ope[i,32] = rand(Bernoulli(reg_grz))
        end

        ##### mortality
        reg_mort = mort_P[sp] * (tanh(6.0 * (mort_reg[sp] - ope[i,5])) + 1)
        ope[33,i] = rand(Bernoulli(reg_mort))

        ##### cell division
        if ope[i,6] ≥ 2 * P_Cquota[sp] * P_Nsuper
            if dvid_type[sp] == 1
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (ope[i,5] - dvid_reg[sp]))+1)
                ope[i,34] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 2
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (ope[i,5] - ope[i,4] - dvid_reg[sp]))+1)
                ope[i,34] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 3
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (ope[i,13] - dvid_reg[sp]))+1)
                ope[i,34] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 4
                cirT = t % 86400 ÷ 3600
                reg_dvid  = dvid_P[sp] * (tanh(dvid_stp[sp] * (cirT - dvid_reg[sp]))+1)
                ope[i,34] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 5
                cirT = t % 86400 ÷ 3600
                regT = dvid_stp[sp]  * (cirT - dvid_reg[sp])
                regS = dvid_stp2[sp] * (ope[i,5] - dvid_reg2[sp])
                reg_dvid  = dvid_P[sp] * (tanh(regS) + 1) * (tanh(regT) + 1)
                ope[i,34] = rand(Bernoulli(reg_dvid))
            else
                throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
            end
        end
    else
        ope[i,32] = 0.0
        ope[i,33] = 0.0
        ope[i,34] = 0.0
    end
end
function probabilities!(ope, sp::AbstractArray{Int64,1}, arch::Architecture, Grz_P, Grz_stp,
                        mort_reg, mort_P, dvid_type, dvid_stp, dvid_stp2, dvid_P,
                        dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    kernel! = probabilities_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, sp, Grz_P, Grz_stp, mort_reg, mort_P, dvid_type,
                    dvid_stp, dvid_stp2, dvid_P, dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ope, inds::AbstractArray{Int64,2},
                                   DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                   lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (ope[i,6] + ope[i,7]) * lossFracC
    POC_con[xi, yi, zi] = POC_con[xi, yi, zi] + (ope[i,6] + ope[i,7]) * (1.0 - lossFracC)
    DON_con[xi, yi, zi] = DON_con[xi, yi, zi] + (ope[i,6] * R_NC + ope[i,8]) * lossFracN
    PON_con[xi, yi, zi] = PON_con[xi, yi, zi] + (ope[i,6] * R_NC + ope[i,8]) * (1.0 - lossFracN)
    DOP_con[xi, yi, zi] = DOP_con[xi, yi, zi] + (ope[i,6] * R_PC + ope[i,9]) * lossFracP
    POP_con[xi, yi, zi] = POP_con[xi, yi, zi] + (ope[i,6] * R_PC + ope[i,9]) * (1.0 - lossFracP)
end
function calc_loss!(ope, inds::AbstractArray{Int64,2}, arch::Architecture,
                    DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_P)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, inds, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ope, inds::AbstractArray{Int64,2},
                                      DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    DIC_con[xi, yi, zi] = DIC_con[xi, yi, zi] + (ope[i,29] - ope[23,i]) * ΔT
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (ope[i,31] - ope[24,i]) * ΔT
    NH4_con[xi, yi, zi] = NH4_con[xi, yi, zi] -  ope[i,25] * ΔT
    NO3_con[xi, yi, zi] = NO3_con[xi, yi, zi] -  ope[i,26] * ΔT
    PO4_con[xi, yi, zi] = PO4_con[xi, yi, zi] -  ope[i,27] * ΔT
end
function calc_consume!(ope, inds::AbstractArray{Int64,2}, arch::Architecture,
                       DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(ope,1),))
    event = kernel!(ope, inds, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g, ΔT)
    wait(device(arch), event)
    return nothing
end
