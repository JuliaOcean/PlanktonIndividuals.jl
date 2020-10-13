##### find and calculate nutrients, αI, and tempfunc for each individual
@kernel function find_NPT_kernel!(plank, inds::AbstractArray{Int64,2},
                                  NH4, NO3, PO4, DOC, par, temp, pop, g::Grids,
                                  α, Φ, TempAe, Tempref, TempCoeff)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds plank[i,16] = max(1.0e-10, NH4[xi+g.Hx, yi+g.Hy, zi+g.Hz])
        @inbounds plank[i,17] = max(1.0e-10, NO3[xi+g.Hx, yi+g.Hy, zi+g.Hz])
        @inbounds plank[i,18] = max(1.0e-10, PO4[xi+g.Hx, yi+g.Hy, zi+g.Hz])
        @inbounds plank[i,19] = max(1.0e-10, DOC[xi+g.Hx, yi+g.Hy, zi+g.Hz])
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
    @inbounds tmp[1:num,1] .= PCmax .* plank[1:num,5] .^ PC_b .* plank[1:num,21] # PCm
    @inbounds tmp[1:num,2] .= exp.(-plank[1:num,20] .* plank[1:num,10] ./ plank[1:num,6] ./ tmp[1:num,1])

    @inbounds plank[1:num,22] .= tmp[1:num,1] .* (1.0 .- tmp[1:num,2]) .* plank[1:num,6] # PS
end

##### calculate DOC uptake rate (mmolC/individual/second)
function calc_VDOC!(plank, tmp, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC, num::Int64)
    @inbounds tmp[1:num,3] .= plank[1:num,7] ./ (plank[1:num,6] .+ plank[1:num,7]) # Qc
    @inbounds tmp[1:num,4] .= max.(0.0, min.(1.0, (Cqmax .- tmp[1:num,3]) ./ (Cqmax - Cqmin))) # inter regulation
    @inbounds tmp[1:num,5] .= plank[1:num,19] ./ (plank[1:num,19] .+ KsatDOC) # KDOC

    @inbounds plank[1:num,23] .= VDOCmax .* plank[1:num,5] .^ VDOC_b .*
        tmp[1:num,4] .* tmp[1:num,5] .* plank[1:num,21] .* plank[1:num,6] # VDOC

    @inbounds plank[1:num,23] .= min.(plank[1:num,19] .* g.V ./ 10.0 ./ ΔT, plank[1:num,23])
end

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
function calc_VN!(plank, tmp, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC, num::Int64)
    @inbounds tmp[1:num,6] .= (plank[1:num,8] .+ plank[1:num,6] .* R_NC) ./ (plank[1:num,6] .+ plank[1:num,7]) # Qn
    @inbounds tmp[1:num,7] .= max.(0.0, min.(1.0, (Nqmax .- tmp[1:num,6]) ./ (Nqmax - Nqmin))) # inter regulation
    @inbounds tmp[1:num,8] .= plank[1:num,16] ./ (plank[1:num,16] .+ KsatNH4) # KNH4
    @inbounds tmp[1:num,9] .= plank[1:num,17] ./ (plank[1:num,17] .+ KsatNO3) # KNO3

    @inbounds plank[1:num,24] .= VNH4max .* plank[1:num,5] .^ VN_b .*
        tmp[1:num,7] .* tmp[1:num,8] .* plank[1:num,21] .* plank[1:num,6] # VNH4

    @inbounds plank[1:num,25] .= VNO3max .* plank[1:num,5] .^ VN_b .*
        tmp[1:num,7] .* tmp[1:num,9] .* plank[1:num,21] .* plank[1:num,6] # VNO3

    @inbounds plank[1:num,24] .= min.(plank[1:num,16] .* g.V ./ 10.0 ./ ΔT, plank[1:num,24])
    @inbounds plank[1:num,25] .= min.(plank[1:num,17] .* g.V ./ 10.0 ./ ΔT, plank[1:num,25])
end

##### calculate PO4 uptake rate (mmolP/individual/second)
function calc_VP!(plank, tmp, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC, num::Int64)
    @inbounds tmp[1:num,10] .= (plank[1:num,9] .+ plank[1:num,6] .* R_PC) ./ (plank[1:num,6] .+ plank[1:num,7]) # Qp
    @inbounds tmp[1:num,11] .= max.(0.0, min.(1.0, (Pqmax .- tmp[1:num,10]) ./ (Pqmax - Pqmin))) # inter regulation
    @inbounds tmp[1:num,12] .= plank[1:num,18] ./ (plank[1:num,18] .+ KsatPO4) # KPO4

    @inbounds plank[1:num,26] .= VPO4max .* plank[1:num,5] .^ VP_b .*
        tmp[1:num,11] .* tmp[1:num,12] .* plank[1:num,21] .* plank[1:num,6] # VPO4

    @inbounds plank[1:num,26] .= min.(plank[1:num,18] .* g.V ./ 10.0 ./ ΔT, plank[1:num,26])
end

##### calculate ρchl
function calc_ρchl!(plank, tmp, Chl2N, num::Int64)
    @inbounds tmp[1:num,13] .= max.(1.0e-10, plank[1:num,20]) .* plank[1:num,10] ./ plank[1:num,6]
    @inbounds tmp[1:num,14] .= isless.(0.1, plank[1:num,20])
    @inbounds plank[1:num,27] .= plank[1:num,22] ./ plank[1:num,6] .* Chl2N ./ tmp[1:num,13]
    @inbounds plank[1:num,27] .= plank[1:num,27] .* tmp[1:num,14]
end

##### calculate respiration (mmolC/individual/second)
function calc_respir!(plank, respir_a, respir_b, num::Int64)
    @inbounds plank[1:num,28] .= respir_a .* plank[1:num,5] .^ respir_b .* plank[1:num,6] .* plank[1:num,21]
end

##### update C, N, P quotas
function update_quotas!(plank, tmp, R_NC, R_PC, ΔT, num::Int64)
    @inbounds plank[1:num,7] .= plank[1:num,7] .+ ΔT .* (plank[1:num,22] .+ plank[1:num,23] .- plank[1:num,28])
    @inbounds plank[1:num,8] .= plank[1:num,8] .+ ΔT .* (plank[1:num,24] .+ plank[1:num,25])
    @inbounds plank[1:num,9] .= plank[1:num,9] .+ ΔT .*  plank[1:num,26]

    @inbounds tmp[1:num,15] .= 0.0 .- plank[1:num,7] # excess
    @inbounds tmp[1:num,15] .= max.(0.0, tmp[1:num,15]) # excess
    @inbounds plank[1:num,7] .= plank[1:num,7] .+ tmp[1:num,15]
    @inbounds plank[1:num,6] .= plank[1:num,6] .- tmp[1:num,15]         # use biomass for respiration
    @inbounds plank[1:num,8] .= plank[1:num,8] .+ tmp[1:num,15] .* R_NC # return N from function pool to N reserve
    @inbounds plank[1:num,9] .= plank[1:num,9] .+ tmp[1:num,15] .* R_PC # return P from function pool to P reserve
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
function calc_BS!(plank, tmp, k_mtb, b_k_mtb, R_NC, R_PC, num::Int64)
    @inbounds tmp[1:num,16] .= k_mtb .* plank[1:num,5] .^ b_k_mtb .* plank[1:num,7]
    @inbounds plank[1:num,29] .= min.(plank[1:num,7], plank[1:num,8] ./ R_NC, plank[1:num,9] ./ R_PC) .*
        k_mtb .* plank[1:num,5] .^ b_k_mtb
    @inbounds plank[1:num,30] .= max.(0.0, tmp[1:num,16] .- plank[1:num,29])
end

##### update C, N, P quotas, biomass, Chla, cell size
function update_biomass!(plank, R_NC, R_PC, Cquota, Nsuper, ΔT, num::Int64)
    @inbounds plank[1:num,6]  .= plank[1:num,6]  .+ ΔT .*  plank[1:num,29]
    @inbounds plank[1:num,7]  .= plank[1:num,7]  .- ΔT .* (plank[1:num,29] .+ plank[1:num,30])
    @inbounds plank[1:num,8]  .= plank[1:num,8]  .- ΔT .*  plank[1:num,29] .* R_NC
    @inbounds plank[1:num,9]  .= plank[1:num,9]  .- ΔT .*  plank[1:num,29] .* R_PC
    @inbounds plank[1:num,10] .= plank[1:num,10] .+ ΔT .*  plank[1:num,29] .* R_NC .* plank[1:num,27]
    @inbounds plank[1:num,12] .= plank[1:num,12] .+ ΔT ./ 3600
    @inbounds plank[1:num,5]  .=(plank[1:num,6]  .+ plank[1:num,7]) ./ Cquota ./ Nsuper
end

##### calculate probability of cell division
##### sizer
function calc_dvid_size!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, num::Int64)
    @inbounds plank[1:num,33] .= dvid_P .* (tanh.(dvid_stp .* (plank[1:num,5] .- dvid_reg)) .+ 1)
    @inbounds plank[1:num,33] .= plank[1:num,33] .* isless.(2*Cquota*Nsuper, plank[1:num,6])
end
##### adder
function calc_dvid_add!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, num::Int64)
    @inbounds plank[1:num,33] .= dvid_P .* (tanh.(dvid_stp .* (plank[1:num,5] .- plank[1:num,4] .- dvid_reg)) .+ 1)
    @inbounds plank[1:num,33] .= plank[1:num,33] .* isless.(2*Cquota*Nsuper, plank[1:num,6])
end
##### age
function calc_dvid_age!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, num::Int64)
    @inbounds plank[1:num,33] .= dvid_P .* (tanh.(dvid_stp .* (plank[1:num,12] .-  dvid_reg)) .+ 1)
    @inbounds plank[1:num,33] .= plank[1:num,33] .* isless.(2*Cquota*Nsuper, plank[1:num,6])
end
##### timer
function calc_dvid_time!(plank, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, t, num::Int64)
    @inbounds plank[1:num,33] .= dvid_P * (tanh(dvid_stp * (t % 86400 ÷ 3600 -  dvid_reg)) + 1)
    @inbounds plank[1:num,33] .= plank[1:num,33] .* isless.(2*Cquota*Nsuper, plank[1:num,6])
end
##### timer & sizer
function calc_dvid_ts!(plank, dvid_stp, dvid_stp2, dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t, num::Int64)
    @inbounds plank[1:num,33] .= dvid_P .* (tanh.(dvid_stp2 .* (plank[1:num,5] .- dvid_reg2)) .+ 1) .*
                                  (tanh(dvid_stp   * (t % 86400 ÷ 3600 - dvid_reg)) + 1)
    @inbounds plank[1:num,33] .= plank[1:num,33] .* isless.(2*Cquota*Nsuper, plank[1:num,6])
end

##### calculate the probability of grazing
##### quadratic grazing
function calc_graz_quadratic!(plank, grz_P, num::Int64)
    @inbounds plank[1:num,31] .= plank[1:num,60] ./ grz_P
end
##### linear grazing decrease with depth
function calc_graz_linear!(plank, grz_P, grz_stp, num::Int64)
    @inbounds plank[1:num,31] .= 1.0 ./ grz_P .* max.(0.15, 1 .- abs.(plank[1:num,3]) ./ grz_stp)
end

##### calculate the probability of mortality
function calc_mort!(plank, mort_reg, mort_P, num::Int64)
    @inbounds plank[1:num,32] .= mort_P .* (tanh.(6.0 .* (mort_reg .- plank[1:num,5])) .+ 1)
end

##### generate the random results from probabilities of grazing, mortality and cell division
function get_rands!(plank, rnd, num::Int64)
    @inbounds plank[1:num,31] .= isless.(rnd[1:num,1], plank[1:num,31])
    @inbounds plank[1:num,32] .= isless.(rnd[1:num,2], plank[1:num,32])
    @inbounds plank[1:num,33] .= isless.(rnd[1:num,3], plank[1:num,33])
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(cts, plank, inds::AbstractArray{Int64,2}, g::Grids, ΔT)
    i = @index(Global, Linear)
    gi = @index(Group)
    if plank[i,58] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds cts[xi, yi, zi, gi, 1] = cts[xi, yi, zi, gi, 1] + (plank[i,28] - plank[i,22]) * ΔT
        @inbounds cts[xi, yi, zi, gi, 5] = cts[xi, yi, zi, gi, 5] + (plank[i,30] - plank[i,23]) * ΔT
        @inbounds cts[xi, yi, zi, gi, 2] = cts[xi, yi, zi, gi, 2] -  plank[i,24] * ΔT
        @inbounds cts[xi, yi, zi, gi, 3] = cts[xi, yi, zi, gi, 3] -  plank[i,25] * ΔT
        @inbounds cts[xi, yi, zi, gi, 4] = cts[xi, yi, zi, gi, 4] -  plank[i,26] * ΔT
    end
end
function calc_consume!(cts, plank, inds::AbstractArray{Int64,2}, arch::Architecture, g::Grids, ΔT)
    kernel! = calc_consume_kernel!(device(arch), 1, (size(plank,1),))
    event = kernel!(cts, plank, inds, g, ΔT)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(cts, plank, inds::AbstractArray{Int64,2}, g::Grids,
                                   lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    if plank[i,58] == 1.0
        @inbounds xi = inds[i,1]
        @inbounds yi = inds[i,2]
        @inbounds zi = inds[i,3]
        @inbounds cts[xi, yi, zi, gi, 5] += (plank[i,6] + plank[i,7]) * lossFracC
        @inbounds cts[xi, yi, zi, gi, 8] += (plank[i,6] + plank[i,7]) * (1.0 - lossFracC)
        @inbounds cts[xi, yi, zi, gi, 6] += (plank[i,6] * R_NC + plank[i,8]) * lossFracN
        @inbounds cts[xi, yi, zi, gi, 9] += (plank[i,6] * R_NC + plank[i,8]) * (1.0 - lossFracN)
        @inbounds cts[xi, yi, zi, gi, 7] += (plank[i,6] * R_PC + plank[i,9]) * lossFracP
        @inbounds cts[xi, yi, zi, gi,10] += (plank[i,6] * R_PC + plank[i,9]) * (1.0 - lossFracP)
    end
end
function calc_loss!(cts, plank, inds::AbstractArray{Int64,2}, arch::Architecture, g::Grids,
                    lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    kernel! = calc_loss_kernel!(device(arch), 1, (size(plank,1),))
    event = kernel!(cts, plank, inds, g, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end
