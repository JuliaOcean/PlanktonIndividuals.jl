##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(plank, inds::AbstractArray{Int64,2},
                                  NH4, NO3, PO4, DOC, par, temp, g::Grids,
                                  α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds plank[i,16] = max(0.0, NH4[xi, yi, zi])
        @inbounds plank[i,17] = max(0.0, NO3[xi, yi, zi])
        @inbounds plank[i,18] = max(0.0, PO4[xi, yi, zi])
        @inbounds plank[i,19] = max(0.0, DOC[xi, yi, zi])
        @inbounds plank[i,20] = α * par[xi, yi, zi] * Φ
        @inbounds plank[i,21] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15)
                                                           - 1.0 / Tempref))) * TempCoeff
    end
end
function find_NPT!(plank, inds::AbstractArray{Int64,2}, arch::Architecture,
                   NH4, NO3, PO4, DOC, par, temp, g::Grids,
                   α, Φ, TempAe, Tempref, TempCoeff)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(plank,1),))

    event = kernel!(plank, inds,  NH4, NO3, PO4, DOC, par, temp, g,
                    α, Φ, TempAe, Tempref, TempCoeff)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@kernel function calc_PS_kernel!(plank, PCmax, PC_b)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds PCm = PCmax * plank[i,5]^PC_b * plank[i,21]
        @inbounds plank[i,22] = PCm * (1.0 - exp(-plank[i,20] * plank[i,10] / plank[i,6] / PCm)) * plank[i,6]
    end
end
function calc_PS!(plank, arch::Architecture, PCmax, PC_b)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, PCmax, PC_b)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
@kernel function calc_VDOC_kernel!(plank, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds Qc = plank[i,7] / (plank[i,6] + plank[i,7])
        @inbounds reg= max(0.0, min(1.0, (Cqmax - Qc) / (Cqmax - Cqmin)))
        @inbounds ksat = plank[i,19] / (plank[i,19] + KsatDOC)
        @inbounds plank[i,23] = VDOCmax * plank[i,5]^VDOC_b * ksat * reg * plank[i,21] * plank[i,6]
        @inbounds plank[i,23] = min(plank[i,19] * g.V / 10.0 / ΔT, plank[i,23])
    end
end
function calc_VDOC!(plank, arch::Architecture,
                    g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    kernel! = calc_VDOC_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, g, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    wait(device(arch), event)
    return nothing
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
@kernel function calc_VN_kernel!(plank, g::Grids,
                                 ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds Qn = (plank[i,8] + plank[i,6] * R_NC) / (plank[i,6] + plank[i,7])
        @inbounds reg= max(0.0, min(1.0, (Nqmax - Qn) / (Nqmax - Nqmin)))
        @inbounds ksatNH = plank[i,16] / (plank[i,16] + KsatNH4)
        @inbounds ksatNO = plank[i,17] / (plank[i,17] + KsatNO3)
        @inbounds plank[i,24] = VNH4max * plank[i,5]^VN_b * ksatNH * reg * plank[i,21] * plank[i,6]
        @inbounds plank[i,24] = min(plank[i,16] * g.V / 10.0 / ΔT, plank[i,24])
        @inbounds plank[i,25] = VNO3max * plank[i,5]^VN_b * ksatNO * reg * plank[i,21] * plank[i,6]
        @inbounds plank[i,25] = min(plank[i,17] * g.V / 10.0 / ΔT, plank[i,25])
    end
end
function calc_VN!(plank, arch::Architecture, g::Grids, ΔT,
                  Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    kernel! = calc_VN_kernel!(device(arch), 256, (size(plank,1),))

    event = kernel!(plank, g, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    wait(device(arch), event)
    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
@kernel function calc_VP_kernel!(plank, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds Qp = (plank[i,9] + plank[i,6] * R_PC) / (plank[i,6] + plank[i,7])
        @inbounds reg= max(0.0, min(1.0, (Pqmax - Qp) / (Pqmax - Pqmin)))
        @inbounds ksat = plank[i,18] / (plank[i,18] + KsatPO4)
        @inbounds plank[i,26] = VPO4max * plank[i,5]^VP_b * ksat * reg * plank[i,21] * plank[i,6]
        @inbounds plank[i,26] = min(plank[i,18] * g.V / 10.0 / ΔT, plank[i,26])
    end
end
function calc_VP!(plank, arch::Architecture, g::Grids,
                  ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    kernel! = calc_VP_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, g, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(plank, Chl2N)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds plank[i,27] = plank[i,20] ≤  0 ? 0.0 :
            plank[i,22] / plank[i,6] * Chl2N / (plank[i,20] * plank[i,10] / plank[i,6])
    end
end
function calc_ρchl!(plank, arch::Architecture, Chl2N)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, respir_a, respir_b)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds plank[i,28] = respir_a * plank[i,5]^respir_b * plank[i,6] * plank[i,21]
    end
end
function calc_respir!(plank, arch::Architecture, respir_a, respir_b)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, respir_a, respir_b)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas
@kernel function update_quotas_kernel!(plank, R_NC, R_PC, ΔT)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds plank[i,7] = plank[i,7] + ΔT * (plank[i,22] + plank[i,23] - plank[i,28])
        if plank[i,7] ≤ 0.0  # if C reserve is not enough for respiration
            @inbounds exceed = 0.0 - plank[i,7]
            @inbounds plank[i,7] = 0.0
            @inbounds plank[i,8] = plank[i,8] + exceed * R_NC # return N from function pool to N reserve
            @inbounds plank[i,9] = plank[i,9] + exceed * R_PC # return P from function pool to P reserve
        end
        @inbounds plank[i,8] = plank[i,8] + ΔT * (plank[i,24] + plank[i,25])
        @inbounds plank[i,9] = plank[i,9] + ΔT *  plank[i,26]
    end
end
function update_quotas!(plank, arch::Architecture, R_NC, R_PC, ΔT)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, R_NC, R_PC, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, k_mtb, b_k_mtb, R_NC, R_PC)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds BS_Cmax = k_mtb * plank[i,5]^b_k_mtb * plank[i,7]
        @inbounds plank[i,29] = min(plank[i,7], plank[i,8]/R_NC, plank[i,9]/R_PC) * k_mtb * plank[i,5]^b_k_mtb
        @inbounds plank[i,30] = max(0.0, BS_Cmax - plank[i,29])
    end
end
function calc_BS!(plank, arch::Architecture, k_mtb, b_k_mtb, R_NC, R_PC)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, k_mtb, b_k_mtb, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
@kernel function update_biomass_kernel!(plank, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds plank[i,6] = plank[i,6]  + ΔT *  plank[i,29]
        @inbounds plank[i,7] = plank[i,7]  - ΔT * (plank[i,29] + plank[i,30])
        @inbounds plank[i,8] = plank[i,8]  - ΔT *  plank[i,29] * R_NC
        @inbounds plank[i,9] = plank[i,9]  - ΔT *  plank[i,29] * R_PC
        @inbounds plank[i,10]= plank[i,10] + ΔT *  plank[i,29] * R_NC * plank[i,27]
        @inbounds plank[i,12]= plank[i,12] + ΔT /3600
        @inbounds plank[i,5] =(plank[i,6]  + plank[i,7]) / P_Cquota / P_Nsuper
    end
end
function update_biomass!(plank, arch::Architecture, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, R_NC, R_PC, P_Cquota, P_Nsuper, ΔT)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(plank, inds::AbstractArray{Int64,2},
                                      DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1] + g.Hx
        @inbounds yi = inds[i,2] + g.Hy
        @inbounds zi = inds[i,3] + g.Hz
        @inbounds DIC_con[xi, yi, zi] = DIC_con[xi, yi, zi] + (plank[i,28] - plank[i,22]) * ΔT
        @inbounds DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (plank[i,30] - plank[i,23]) * ΔT
        @inbounds NH4_con[xi, yi, zi] = NH4_con[xi, yi, zi] -  plank[i,24] * ΔT
        @inbounds NO3_con[xi, yi, zi] = NO3_con[xi, yi, zi] -  plank[i,25] * ΔT
        @inbounds PO4_con[xi, yi, zi] = PO4_con[xi, yi, zi] -  plank[i,26] * ΔT
    end
end
function calc_consume!(plank, inds::AbstractArray{Int64,2}, arch::Architecture,
                       DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, inds, DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate probability of cell division
@kernel function calc_dvid_kernel!(plank, dvid_type, dvid_stp, dvid_stp2,
                                   dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            if dvid_type == 1
                @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,5] - dvid_reg))+1)
            elseif dvid_type == 2
                @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,5] - plank[i,4] - dvid_reg))+1)
            elseif dvid_type == 3
                @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,12] - dvid_reg))+1)
            elseif dvid_type == 4
                cirT = t % 86400 ÷ 3600
                @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (cirT - dvid_reg))+1)
            elseif dvid_type == 5
                cirT = t % 86400 ÷ 3600
                regT = dvid_stp  * (cirT - dvid_reg)
                @inbounds regS = dvid_stp2 * (plank[i,5] - dvid_reg2)
                @inbounds plank[i,33] = dvid_P * (tanh(regS) + 1) * (tanh(regT) + 1)
            else
                throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
            end
        end
    end
end
function calc_dvid!(plank, arch::Architecture, dvid_type, dvid_stp, dvid_stp2,
                    dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    kernel! = calc_dvid_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_type, dvid_stp, dvid_stp2, dvid_P,
                    dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    wait(device(arch), event)
    return nothing
end
##### generate the random results from probabilities of grazing, mortality and cell division
@kernel function get_rands_kernel!(plank)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds plank[i,31] = plank[i,58] ≥ plank[i,31] ? 0.0 : 1.0
        @inbounds plank[i,32] = plank[i,59] ≥ plank[i,32] ? 0.0 : 1.0
        @inbounds plank[i,33] = plank[i,60] ≥ plank[i,33] ? 0.0 : 1.0
    end
end
function get_rands!(plank, arch::Architecture)
    kernel! = get_rands_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(plank, inds::AbstractArray{Int64,2},
                                   DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                   lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    if plank[i,61] == 1.0
        @inbounds xi = inds[i,1] + g.Hx
        @inbounds yi = inds[i,2] + g.Hy
        @inbounds zi = inds[i,3] + g.Hz
        @inbounds DOC_con[xi, yi, zi] = DOC_con[xi, yi, zi] + (plank[i,6] + plank[i,7]) * lossFracC
        @inbounds POC_con[xi, yi, zi] = POC_con[xi, yi, zi] + (plank[i,6] + plank[i,7]) * (1.0 - lossFracC)
        @inbounds DON_con[xi, yi, zi] = DON_con[xi, yi, zi] + (plank[i,6] * R_NC + plank[i,8]) * lossFracN
        @inbounds PON_con[xi, yi, zi] = PON_con[xi, yi, zi] + (plank[i,6] * R_NC + plank[i,8]) * (1.0 - lossFracN)
        @inbounds DOP_con[xi, yi, zi] = DOP_con[xi, yi, zi] + (plank[i,6] * R_PC + plank[i,9]) * lossFracP
        @inbounds POP_con[xi, yi, zi] = POP_con[xi, yi, zi] + (plank[i,6] * R_PC + plank[i,9]) * (1.0 - lossFracP)
    end
end
function calc_loss!(plank, inds::AbstractArray{Int64,2}, arch::Architecture,
                    DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, inds, DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end
