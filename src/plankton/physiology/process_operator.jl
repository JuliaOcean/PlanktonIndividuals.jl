##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(plank, inds::AbstractArray{Int64,2},
                                  NH4, NO3, PO4, DOC, par, temp, pop, g::Grids,
                                  α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds plank[i,16] = max(1.0e-10, NH4[xi, yi, zi])
        @inbounds plank[i,17] = max(1.0e-10, NO3[xi, yi, zi])
        @inbounds plank[i,18] = max(1.0e-10, PO4[xi, yi, zi])
        @inbounds plank[i,19] = max(1.0e-10, DOC[xi, yi, zi])
        @inbounds plank[i,20] = α * par[xi, yi, zi] * Φ
        @inbounds plank[i,21] = max(1.0e-10, exp(TempAe * (1.0 / (temp[xi, yi, zi] + 273.15)
                                                           - 1.0 / Tempref))) * TempCoeff
        @inbounds plank[i,60] = pop[xi, yi, zi]
    end
end
function find_NPT!(plank, inds::AbstractArray{Int64,2}, arch::Architecture,
                   NH4, NO3, PO4, DOC, par, temp, pop, g::Grids,
                   α, Φ, TempAe, Tempref, TempCoeff)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(plank,1),))

    event = kernel!(plank, inds,  NH4, NO3, PO4, DOC, par, temp, pop, g,
                    α, Φ, TempAe, Tempref, TempCoeff)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
function calc_PS!(plank, tmp, PCmax, PC_b, num::Int64)
    tmp[1:num,1] .= PCmax .* plank[1:num,5] .^ PC_b .* plank[1:num,21] # PCm
    tmp[1:num,2] .= exp.(-plank[1:num,20] .* plank[1:num,10] ./ plank[1:num,6] ./ tmp[1:num,1]) # exp(-αI⋅Chl/C.PCm)

    plank[1:num,22] .= tmp[1:num,1] .* (1.0 .- tmp[1:num,2]) * plank[1:num,6] # PS
end

##### calculate DOC uptake rate (mmolC/individual/second)
function calc_VDOC!(plank, tmp, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    tmp[1:num,3] .= plank[1:num,7] ./ (plank[1:num,6] .+ plank[1:num,7]) # Qc
    tmp[1:num,4] .= max.(0.0, min.(1.0, (Cqmax .- tmp[1:num,3]) ./ (Cqmax - Cqmin))) # intracellular regulation
    tmp[1:num,5] .= plank[1:num,19] ./ (plank[1:num,19] .+ KsatDOC) # KDOC

    plank[1:num,23] .= VDOCmax .* plank[1:num,5] .^ VDOC_b .*
        tmp[1:num,4] .* tmp[1:num,5] * plank[1:num,21] * plank[1:num,6] # VDOC

    plank[1:num,23] .= min.(plank[1:num,19] .* g.V ./ 10.0 ./ ΔT, plank[1:num,23])
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
function calc_VN!(plank, tmp, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    tmp[1:num,6] .= (plank[1:num,8] .+ plank[1:num,6] .* R_NC) ./ (plank[1:num,6] .+ plank[1:num,7]) # Qn
    tmp[1:num,7] .= max.(0.0, min.(1.0, (Nqmax .- tmp[1:num,6]) ./ (Cqmax - Cqmin))) # intracellular regulation
    tmp[1:num,8] .= plank[1:num,16] ./ (plank[1:num,16] .+ KsatNH4) # KNH4
    tmp[1:num,9] .= plank[1:num,17] ./ (plank[1:num,17] .+ KsatNO3) # KNO3

    plank[1:num,24] .= VNH4max .* plank[1:num,5] .^ VN_b .*
        tmp[1:num,7] .* tmp[1:num,8] * plank[1:num,21] * plank[1:num,6] # VNH4

    plank[1:num,25] .= VNO3max .* plank[1:num,5] .^ VN_b .*
        tmp[1:num,7] .* tmp[1:num,9] * plank[1:num,21] * plank[1:num,6] # VNO3

    plank[1:num,24] .= min.(plank[1:num,16] .* g.V ./ 10.0 ./ ΔT, plank[1:num,24])
    plank[1:num,25] .= min.(plank[1:num,17] .* g.V ./ 10.0 ./ ΔT, plank[1:num,25])
end

##### calculate PO4 uptake rate (mmolP/individual/second)
function calc_VP!(plank, tmp, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    tmp[1:num,10] .= (plank[1:num,9] .+ plank[1:num,6] .* R_PC) ./ (plank[1:num,6] .+ plank[1:num,7]) # Qp
    tmp[1:num,11] .= max.(0.0, min.(1.0, (Pqmax .- tmp[1:num,10]) ./ (Pqmax - Pqmin))) # intracellular regulation
    tmp[1:num,12] .= plank[1:num,18] ./ (plank[1:num,18] .+ KsatPO4) # KPO4

    plank[1:num,26] .= VPO4max .* plank[1:num,5] .^ VP_b .*
        tmp[1:num,11] .* tmp[1:num,12] * plank[1:num,21] * plank[1:num,6] # VPO4

    plank[1:num,26] .= min.(plank[1:num,18] .* g.V ./ 10.0 ./ ΔT, plank[1:num,26])
end

##### calculate ρchl
function calc_ρchl!(plank, tmp, Chl2N)
    tmp[1:num,13] .= max.(1.0e-20, plank[1:num,20]) .* plank[1:num,10] ./ plank[1:num,6]
    tmp[1:num,14] .= isless.(0.1, plank[1:num,20])
    plank[1:num,27] .= plank[1:num,22] ./ plank[1:num,6] .* Chl2N ./ tmp[1:num,13]
    plank[1:num,27] .= plank[1:num,27] .* tmp[1:num,14]
end

##### calculate respiration (mmolC/individual/second)
function calc_respir!(plank, tmp, respir_a, respir_b)
    plank[1:num,28] .= respir_a .* plank[1:num,5] .^ respir_b .* plank[1:num,6] .* plank[1:num,21]
end

##### update C, N, P quotas
function update_quotas!(plank, tmp, R_NC, R_PC, ΔT)
    plank[1:num,7] .= plank[1:num,7] .+ ΔT .* (plank[1:num,22] .+ plank[1:num,23] .- plank[1:num,28])
    plank[1:num,8] .= plank[1:num,8] .+ ΔT .* (plank[1:num,24] .+ plank[1:num,25])
    plank[1:num,9] .= plank[1:num,9] .+ ΔT .*  plank[1:num,26]

    tmp[1:num,15] .= 0.0 .- plank[1:num,7] # excess
    tmp[1:num,15] .= max.(0.0, tmp[1:num,15]) # excess
    plank[1:num,7] .= plank[1:num,7] .+ tmp[1:num,15]
    plank[1:num,6] .= plank[1:num,6] .- tmp[1:num,15]         # use biomass for respiration
    plank[1:num,8] .= plank[1:num,8] .+ tmp[1:num,15] .* R_NC # return N from function pool to N reserve
    plank[1:num,9] .= plank[1:num,9] .+ tmp[1:num,15] .* R_PC # return P from function pool to P reserve
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
function calc_BS!(plank, tmp, k_mtb, b_k_mtb, R_NC, R_PC)
    tmp[1:num,16] .= k_mtb .* plank[1:num,5] .^ b_k_mtb .* plank[1:num,7]
    plank[1:num,29] .= min.(plank[1:num,7], plank[1:num,8] ./ R_NC, plank[1:num,9] ./ R_PC) .*
        k_mtb .* plank[1:num,5] .^ b_k_mtb
    plank[1:num,30] .= max.(0.0, tmp[1:num,16] .- plank[1:num,29])
end

##### update C, N, P quotas, biomass, Chla, cell size
function update_biomass!(plank, tmp, R_NC, R_PC, Cquota, Nsuper, ΔT)
    plank[1:num,6]  .= plank[1:num,6]  .+ ΔT .*  plank[1:num,29]
    plank[1:num,7]  .= plank[1:num,7]  .- ΔT .* (plank[1:num,29] .+ plank[1:num,30])
    plank[1:num,8]  .= plank[1:num,8]  .- ΔT .*  plank[1:num,29] .* R_NC
    plank[1:num,9]  .= plank[1:num,9]  .- ΔT .*  plank[1:num,29] .* R_PC
    plank[1:num,10] .= plank[1:num,10] .+ ΔT .*  plank[1:num,29] .* R_NC .* plank[1:num,27]
    plank[1:num,12] .= plank[1:num,12] .+ ΔT ./ 3600
    plank[1:num,5]  .=(plank[1:num,6]  .+ plank[1:num,7]) ./ Cquota ./ Nsuper
end

##### calculate probability of cell division
##### sizer
@kernel function calc_dvid_size_kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,5] - dvid_reg))+1)
        end
    end
end
function calc_dvid_size!(plank, arch::Architecture, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    kernel! = calc_dvid_size_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    wait(device(arch), event)
    return nothing
end
##### adder
@kernel function calc_dvid_add_kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,5] - plank[i,4] - dvid_reg))+1)
        end
    end
end
function calc_dvid_add!(plank, arch::Architecture, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    kernel! = calc_dvid_add_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    wait(device(arch), event)
    return nothing
end
##### age
@kernel function calc_dvid_age_kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (plank[i,12] - dvid_reg))+1)
        end
    end
end
function calc_dvid_age!(plank, arch::Architecture, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    kernel! = calc_dvid_age_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    wait(device(arch), event)
    return nothing
end
##### timer
@kernel function calc_dvid_time_kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, t)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp * (t % 86400 ÷ 3600 - dvid_reg))+1)
        end
    end
end
function calc_dvid_time!(plank, arch::Architecture, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, t)
    kernel! = calc_dvid_time_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, t)
    wait(device(arch), event)
    return nothing
end
##### timer & sizer
@kernel function calc_dvid_ts_kernel!(plank, dvid_stp, dvid_stp2,
                                      dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        if plank[i,6] ≥ 2 * Cquota * Nsuper
            @inbounds plank[i,33] = dvid_P * (tanh(dvid_stp2 * (plank[i,5] - dvid_reg2)) + 1) *
                                             (tanh(dvid_stp  * (t % 86400 ÷ 3600 - dvid_reg)) + 1)
        end
    end
end
function calc_dvid_ts!(plank, arch::Architecture, dvid_stp, dvid_stp2, dvid_P,
                       dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    kernel! = calc_dvid_ts_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, dvid_stp, dvid-stp2, dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    wait(device(arch), event)
    return nothing
end

##### calculate the probability of grazing
##### quadratic grazing
@kernel function calc_graz_quadratic_kernel!(plank, grz_P)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds plank[i,31] = plank[i,60] /grz_P
    end
end
function calc_graz_quadratic!(plank, arch::Architecture, grz_P)
    kernel! = calc_graz_quadratic_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, grz_P)
    wait(device(arch), event)
    return nothing
end
##### linear grazing decrease with depth
@kernel function calc_graz_linear_kernel!(plank, grz_P, grz_stp)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds plank[i,31] = 1.0 / grz_P * max(0.15, 1 - abs(plank[i,3]) / grz_stp)
    end
end
function calc_graz_linear!(plank, arch::Architecture, grz_P, grz_stp)
    kernel! = calc_graz_linear_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, grz_P, grz_stp)
    wait(device(arch), event)
    return nothing
end

##### calculate the probability of mortality
@kernel function calc_mort_kernel!(plank, mort_reg, mort_P)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds plank[i,32] = mort_P * (tanh(6.0 * (mort_reg - plank[i,5]))+1)
    end
end
function calc_mort!(plank, arch::Architecture, mort_reg, mort_P)
    kernel! = calc_mort_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, mort_reg, mort_P)
    wait(device(arch), event)
    return nothing
end

##### generate the random results from probabilities of grazing, mortality and cell division
@kernel function get_rands_kernel!(plank, rnd)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds plank[i,31] = rnd[i,1] ≥ plank[i,31] ? 0.0 : 1.0
        @inbounds plank[i,32] = rnd[i,2] ≥ plank[i,32] ? 0.0 : 1.0
        @inbounds plank[i,33] = rnd[i,3] ≥ plank[i,33] ? 0.0 : 1.0
    end
end
function get_rands!(plank, rnd, arch::Architecture)
    kernel! = get_rands_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, rnd)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(plank, inds::AbstractArray{Int64,2},
                                      DIC_con, DOC_con, NH4_con, NO3_con, PO4_con, g::Grids, ΔT)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
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

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(plank, inds::AbstractArray{Int64,2},
                                   DOC_con, POC_con, DON_con, PON_con, DOP_con, POP_con, g::Grids,
                                   lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
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
