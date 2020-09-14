##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(plank, sp::Int64, inds::AbstractArray{Int64,2},
                                  NH4, NO3, PO4, DOC, par, temp, g::Grids,
                                  α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    plank[i,16] = max(0.0, NH4[xi, yi, zi])
    plank[i,17] = max(0.0, NO3[xi, yi, zi])
    plank[i,18] = max(0.0, PO4[xi, yi, zi])
    plank[i,19] = max(0.0, DOC[xi, yi, zi])
    plank[i,20] = α[sp] * par[xi, yi, zi] * Φ[sp]
    plank[i,21] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15) - 1.0 / Tempref))) * TempCoeff
end
function find_NPT!(plank, sp::Int64, inds::AbstractArray{Int64,2}, arch::Architecture,
                   NH4, NO3, PO4, DOC, par, temp, g::Grids,
                   α, Φ, TempAe, Tempref, TempCoeff)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(plank,1),))

    event = kernel!(plank, sp, inds,  NH4, NO3, PO4, DOC, par, temp, g,
                    α, Φ, TempAe, Tempref, TempCoeff)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@kernel function calc_PS_kernel!(plank, sp::Int64, PCmax, PC_b)
    i = @index(Global, Linear)
    PCm = PCmax[sp] * plank[i,5]^PC_b[sp] * plank[i,21]
    plank[i,22] = PCm * (1.0 - exp(-plank[i,20] * plank[i,10] / plank[i,6] / PCm)) * plank[i,6]
end
function calc_PS!(plank, sp::Int64, arch::Architecture, PCmax, PC_b)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, PCmax, PC_b)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
@kernel function calc_VDOC_kernel!(plank, sp::Int64, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    i = @index(Global, Linear)

    Qc = plank[i,7] / (plank[i,6] + plank[i,7])
    reg= max(0.0, min(1.0, (Cqmax[sp] - Qc) / (Cqmax[sp] - Cqmin[sp])))

    ksat = plank[i,19] / (plank[i,19] + KsatDOC[sp])

    plank[i,23] = VDOCmax[sp] * plank[i,5]^VDOC_b[sp] * ksat * reg * plank[i,21] * plank[i,6]
    plank[i,23] = min(plank[i,19] * g.V / 10.0 / ΔT, plank[i,23])
end
function calc_VDOC!(plank, sp::Int64, arch::Architecture,
                    g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    kernel! = calc_VDOC_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, g, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    wait(device(arch), event)
    return nothing
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
@kernel function calc_VN_kernel!(plank, sp::Int64, g::Grids,
                                 ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    i = @index(Global, Linear)

    Qn = (plank[i,8] + plank[i,6] * R_NC[sp]) / (plank[i,6] + plank[i,7])
    reg= max(0.0, min(1.0, (Nqmax[sp] - Qn) / (Nqmax[sp] - Nqmin[sp])))

    ksatNH = plank[i,16] / (plank[i,16] + KsatNH4[sp])
    ksatNO = plank[i,17] / (plank[i,17] + KsatNO3[sp])

    plank[i,24] = VNH4max[sp] * plank[i,5]^VN_b[sp] * ksatNH * reg * plank[i,21] * plank[i,6]
    plank[i,24] = min(plank[i,16] * g.V / 10.0 / ΔT, plank[i,24])

    plank[i,25] = VNO3max[sp] * plank[i,5]^VN_b[sp] * ksatNO * reg * plank[i,21] * plank[i,6]
    plank[i,25] = min(plank[i,17] * g.V / 10.0 / ΔT, plank[i,25])
end
function calc_VN!(plank, sp::Int64, arch::Architecture, g::Grids, ΔT,
                  Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    kernel! = calc_VN_kernel!(device(arch), 256, (size(plank,1),))

    event = kernel!(plank, sp, g, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    wait(device(arch), event)
    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
@kernel function calc_VP_kernel!(plank, sp::Int64, g::Grids,
                                 ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    i = @index(Global, Linear)

    Qp = (plank[i,9] + plank[i,6] * R_PC[sp]) / (plank[i,6] + plank[i,7])
    reg= max(0.0, min(1.0, (Pqmax[sp] - Qp) / (Pqmax[sp] - Pqmin[sp])))

    ksat = plank[i,18] / (plank[i,18] + KsatPO4[sp])

    plank[i,26] = VPO4max[sp] * plank[i,5]^VP_b[sp] * ksat * reg * plank[i,21] * plank[i,6]
    plank[i,26] = min(plank[i,18] * g.V / 10.0 / ΔT, plank[i,26])
end
function calc_VP!(plank, sp::Int64, arch::Architecture, g::Grids,
                  ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    kernel! = calc_VP_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, g, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(plank, sp::Int64, Chl2N)
    i = @index(Global, Linear)
    plank[i,27] = plank[i,20] < 0 ? 0.0 :
        plank[i,22] / plank[i,6] * Chl2N[sp] / (plank[i,20] * plank[i,10] / plank[i,6])
end
function calc_ρchl!(plank, sp::Int64, arch::Architecture, Chl2N)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, sp::Int64, respir_a, respir_b)
    i = @index(Global, Linear)
    plank[i,28] = respir_a[sp] * plank[i,5]^respir_b[sp] * plank[i,6] * plank[i,21]
end
function calc_respir!(plank, sp::Int64, arch::Architecture, respir_a, respir_b)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, respir_a, respir_b)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas
@kernel function update_quotas_kernel!(plank, sp::Int64, R_NC, R_PC, ΔT)
    i = @index(Global, Linear)
    plank[i,7] = plank[i,7] + ΔT * (plank[i,22] + plank[i,23] - plank[i,28])
    if plank[i,7] ≤ 0.0  # if C reserve is not enough for respiration
        exceed = 0.0 - plank[i,7]
        plank[i,7] = 0.0
        plank[i,8] = plank[i,8] + exceed * R_NC[sp] # return N from function pool to N reserve
        plank[i,9] = plank[i,9] + exceed * R_PC[sp] # return P from function pool to P reserve
    end
    plank[i,8] = plank[i,8] + ΔT * (plank[i,24] + plank[i,25])
    plank[i,9] = plank[i,9] + ΔT *  plank[i,26]
end
function update_quotas!(plank, sp::Int64, arch::Architecture, R_NC, R_PC, ΔT)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, R_NC, R_PC, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, sp::Int64, k_mtb, b_k_mtb, R_NC, R_PC)
    i = @index(Global, Linear)
    BS_Cmax = k_mtb[sp] * plank[i,5]^b_k_mtb[sp] * plank[i,7]
    plank[i,29] = min(plank[i,7], plank[i,8]/R_NC[sp], plank[i,9]/R_PC[sp]) * k_mtb[sp] * plank[i,5]^b_k_mtb[sp]
    plank[i,30] = max(0.0, BS_Cmax - plank[i,29])
end
function calc_BS!(plank, sp::Int64, arch::Architecture, k_mtb, b_k_mtb, R_NC, R_PC)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, k_mtb, b_k_mtb, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
@kernel function update_biomass_kernel!(plank, sp::Int64, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    i = @index(Global, Linear)
    plank[i,6] = plank[i,6]  + ΔT *  plank[i,29]
    plank[i,7] = plank[i,7]  - ΔT * (plank[i,29] + plank[i,30])
    plank[i,8] = plank[i,8]  - ΔT *  plank[i,29] * R_NC[sp]
    plank[i,9] = plank[i,9]  - ΔT *  plank[i,29] * R_PC[sp]
    plank[i,10]= plank[i,10] + ΔT *  plank[i,29] * R_NC[sp] * plank[i,27]
    plank[i,12]= plank[i,12] + ΔT /3600
    plank[i,5] =(plank[i,6]  + plank[i,7]) / P_Cquota[sp] / P_Nsuper[sp]

end
function update_biomass!(plank, sp::Int64, arch::Architecture, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of grazing, mortality and cell division
@kernel function probabilities_kernel!(plank, sp::Int64, Grz_P, Grz_stp,
                                       mort_reg, mort_P, dvid_type, dvid_stp, dvid_stp2,
                                       dvid_P, dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    i = @index(Global, Linear)

    if t%300 == 1
        ##### grazing
        if Grz_P[sp] == 0
            plank[i,31] = 0.0
        else
            reg_grz = Grz_stp[sp] == 0 ? 1.0 / Grz_P[sp] :
                1.0 / Grz_P[sp] * max(0.15, 1 - abs(plank[i,3]) / Grz_stp[sp])
            # reg = 1.0 / Grz_P
            plank[i,31] = rand(Bernoulli(reg_grz))
        end

        ##### mortality
        reg_mort = mort_P[sp] * (tanh(6.0 * (mort_reg[sp] - plank[i,5])) + 1)
        plank[32,i] = rand(Bernoulli(reg_mort))

        ##### cell division
        if plank[i,6] ≥ 2 * P_Cquota[sp] * P_Nsuper
            if dvid_type[sp] == 1
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (plank[i,5] - dvid_reg[sp]))+1)
                plank[i,33] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 2
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (plank[i,5] - plank[i,4] - dvid_reg[sp]))+1)
                plank[i,33] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 3
                reg_dvid = dvid_P[sp] * (tanh(dvid_stp[sp] * (plank[i,12] - dvid_reg[sp]))+1)
                plank[i,33] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 4
                cirT = t % 86400 ÷ 3600
                reg_dvid  = dvid_P[sp] * (tanh(dvid_stp[sp] * (cirT - dvid_reg[sp]))+1)
                plank[i,33] = Float64(rand(Bernoulli(reg_dvid)))
            elseif dvid_type[sp] == 5
                cirT = t % 86400 ÷ 3600
                regT = dvid_stp[sp]  * (cirT - dvid_reg[sp])
                regS = dvid_stp2[sp] * (plank[i,5] - dvid_reg2[sp])
                reg_dvid  = dvid_P[sp] * (tanh(regS) + 1) * (tanh(regT) + 1)
                plank[i,33] = rand(Bernoulli(reg_dvid))
            else
                throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
            end
        end
    else
        plank[i,31] = 0.0
        plank[i,32] = 0.0
        plank[i,33] = 0.0
    end
end
function probabilities!(plank, sp::Int64, arch::Architecture, Grz_P, Grz_stp,
                        mort_reg, mort_P, dvid_type, dvid_stp, dvid_stp2, dvid_P,
                        dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    kernel! = probabilities_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, Grz_P, Grz_stp, mort_reg, mort_P, dvid_type,
                    dvid_stp, dvid_stp2, dvid_P, dvid_reg, dvid_reg2, P_Cquota, P_Nsuper, t)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(plank, sp::Int64, inds::AbstractArray{Int64,2},
                                   DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                   lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (plank[i,6] + plank[i,7]) * lossFracC[sp]
    POC_con[xi, yi, zi] = POC_con[xi, yi, zi] + (plank[i,6] + plank[i,7]) * (1.0 - lossFracC[sp])
    DON_con[xi, yi, zi] = DON_con[xi, yi, zi] + (plank[i,6] * R_NC[sp] + plank[i,8]) * lossFracN[sp]
    PON_con[xi, yi, zi] = PON_con[xi, yi, zi] + (plank[i,6] * R_NC[sp] + plank[i,8]) * (1.0 - lossFracN[sp])
    DOP_con[xi, yi, zi] = DOP_con[xi, yi, zi] + (plank[i,6] * R_PC[sp] + plank[i,9]) * lossFracP[sp]
    POP_con[xi, yi, zi] = POP_con[xi, yi, zi] + (plank[i,6] * R_PC[sp] + plank[i,9]) * (1.0 - lossFracP[sp])
end
function calc_loss!(plank, sp::Int64, inds::AbstractArray{Int64,2}, arch::Architecture,
                    DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_P)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, sp, inds, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(plank, inds::AbstractArray{Int64,2},
                                      DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    DIC_con[xi, yi, zi] = DIC_con[xi, yi, zi] + (plank[i,28] - plank[22,i]) * ΔT
    DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (plank[i,30] - plank[23,i]) * ΔT
    NH4_con[xi, yi, zi] = NH4_con[xi, yi, zi] -  plank[i,24] * ΔT
    NO3_con[xi, yi, zi] = NO3_con[xi, yi, zi] -  plank[i,25] * ΔT
    PO4_con[xi, yi, zi] = PO4_con[xi, yi, zi] -  plank[i,26] * ΔT
end
function calc_consume!(plank, inds::AbstractArray{Int64,2}, arch::Architecture,
                       DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, inds, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g, ΔT)
    wait(device(arch), event)
    return nothing
end
