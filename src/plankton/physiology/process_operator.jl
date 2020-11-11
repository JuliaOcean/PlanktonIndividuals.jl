##### find and calculate nutrients, αI, and tempfunc for each individual
function find_NPT!(nuts, x, y, z, ac, g::Grids, NH4, NO3, PO4, DOC, par, temp, pop, α, Φ, TempAe, Tempref, TempCoeff)
    @inbounds nuts.NH4 .= NH4[CartesianIndex.(x .+ g.Hx, y .+ g.Hy, z .+ g.Hz)] .* ac
    @inbounds nuts.NO3 .= NO3[CartesianIndex.(x .+ g.Hx, y .+ g.Hy, z .+ g.Hz)] .* ac
    @inbounds nuts.PO4 .= PO4[CartesianIndex.(x .+ g.Hx, y .+ g.Hy, z .+ g.Hz)] .* ac
    @inbounds nuts.DOC .= DOC[CartesianIndex.(x .+ g.Hx, y .+ g.Hy, z .+ g.Hz)] .* ac
    @inbounds nuts.αI  .= par[CartesianIndex.(x, y, z)] .* α .* Φ .* ac
    @inbounds nuts.Tem .= temp[CartesianIndex.(x, y, z)] .* ac
    @inbounds nuts.pop .= pop[CartesianIndex.(x, y, z)] .* ac

    @inbounds nuts.NH4 .= max.(1.0e-10, nuts.NH4) .* ac
    @inbounds nuts.NO3 .= max.(1.0e-10, nuts.NO3) .* ac
    @inbounds nuts.PO4 .= max.(1.0e-10, nuts.PO4) .* ac
    @inbounds nuts.DOC .= max.(1.0e-10, nuts.DOC) .* ac
    @inbounds nuts.Tem .= max.(1.0e-10, exp.(TempAe .* (1.0 ./ (nuts.Tem .+ 273.15) .-
                                                        (1.0/Tempref)))) .* TempCoeff .* ac

    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
function calc_PS!(plank, proc, nuts, PCmax, PC_b)
    @inbounds proc.PS .= PCmax .* plank.Sz .^ PC_b .* nuts.Tem .* plank.ac
    @inbounds proc.PS .= proc.PS .* (1.0 .- exp.(-nuts.αI .* plank.chl ./
                                                 max.(1.0e-10, plank.Bm .* proc.PS))) .* plank.Bm
    @inbounds proc.PS .= proc.PS .* plank.ac

    return nothing
end

##### calculate the intracellular regulations of nutrient uptake
@inline regQ(Nq, Bm, Cq, Nqmax, Nqmin, R_NC) = max(0.0, min(1.0, (Nqmax - (Nq + Bm * R_NC) /
                                                                  max(1.0e-10, Bm + Cq)) / (Nqmax - Nqmin)))

##### calculate NH4 and NO3 uptake rate (mmolN/individual/second)
function calc_VN!(plank, proc, nuts, g::Grids, ΔT, Nqmax, Nqmin, VNH4max, VNO3max, VN_b, KsatNH4, KsatNO3, R_NC)
    @inbounds proc.VNH4 .= VNH4max .* plank.Sz .^ VN_b .*
        regQ.(plank.Nq, plank.Bm, plank.Cq, Nqmax, Nqmin, R_NC) .*
        nuts.NH4 ./ (nuts.NH4 .+ KsatNH4) .* nuts.Tem .* plank.Bm
    @inbounds proc.VNH4 .= min.(nuts.NH4 .* g.V ./10.0 ./ ΔT, proc.VNH4) .* plank.ac

    @inbounds proc.VNO3 .= VNO3max .* plank.Sz .^ VN_b .*
        regQ.(plank.Nq, plank.Bm, plank.Cq, Nqmax, Nqmin, R_NC) .*
        nuts.NO3 ./ (nuts.NO3 .+ KsatNO3) .* nuts.Tem .* plank.Bm
    @inbounds proc.VNO3 .= min.(nuts.NO3 .* g.V ./10.0 ./ ΔT, proc.VNO3) .* plank.ac

    return nothing
end

##### calculate PO4 uptake rate (mmolP/individual/second)
function calc_VP!(plank, proc, nuts, g::Grids, ΔT, Pqmax, Pqmin, VPO4max, VP_b, KsatPO4, R_PC)
    @inbounds proc.VPO4 .= VPO4max .* plank.Sz .^ VP_b .*
        regQ.(plank.Pq, plank.Bm, plank.Cq, Pqmax, Pqmin, R_PC) .*
        nuts.PO4 ./ (nuts.PO4 .+ KsatPO4) .* nuts.Tem .* plank.Bm
    @inbounds proc.VPO4 .= min.(nuts.PO4 .* g.V ./10.0 ./ ΔT, proc.VPO4) .* plank.ac

    return nothing
end

##### calculate ρchl
function calc_ρchl!(plank, proc, nuts, Chl2N)
    @inbounds proc.ρchl .= proc.PS ./ max.(1.0e-10, plank.Bm) .*
        Chl2N ./ max.(1.0e-10, nuts.αI .* plank.chl ./ plank.Bm)
    @inbounds proc.ρchl .= proc.ρchl .* isless.(0.1, nuts.αI) .* plank.ac

    return nothing
end

##### calculate respiration (mmolC/individual/second)
function calc_respir!(plank, proc, nuts, respir_a, respir_b)
    @inbounds proc.resp .= respir_a .* plank.Sz .^ respir_b .* plank.Bm .* nuts.Tem .* plank.ac

    return nothing
end

##### update C, N, P quotas for the first time of each time step
function update1_quotas!(plank, proc, ΔT)
    @inbounds plank.Cq .= plank.Cq .+ ΔT .*  proc.PS
    @inbounds plank.Nq .= plank.Nq .+ ΔT .* (proc.VNH4 .+ proc.VNO3)
    @inbounds plank.Pq .= plank.Pq .+ ΔT .* proc.VPO4

    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
function calc_VDOC!(plank, proc, nuts, g::Grids, ΔT, Cqmax, Cqmin, VDOCmax, VDOC_b, KsatDOC)
    @inbounds proc.VDOC .= VDOCmax .* plank.Sz .^ VDOC_b .*
        regQ.(plank.Cq, plank.Bm, plank.Cq, Cqmax, Cqmin, 0.0) .*
        nuts.DOC ./ (nuts.DOC .+ KsatDOC) .* nuts.Tem .* plank.Bm
    @inbounds proc.VDOC .= min.(nuts.DOC .* g.V ./10.0 ./ ΔT, proc.VDOC) .* plank.ac

    return nothing
end

##### update C, N, P quotas for the second time of each time step
function update2_quotas!(plank, proc, ΔT, R_NC, R_PC)
    @inbounds plank.Cq .= plank.Cq .+ ΔT .* (proc.VDOC .- proc.resp)
    @inbounds plank.Cq .= plank.Cq .+ max.(0.0, (0.0 .- plank.Cq))
    @inbounds plank.Bm .= plank.Bm .- max.(0.0, (0.0 .- plank.Cq))
    @inbounds plank.Nq .= plank.Nq .+ max.(0.0, (0.0 .- plank.Cq)) .* R_NC
    @inbounds plank.Pq .= plank.Pq .+ max.(0.0, (0.0 .- plank.Cq)) .* R_PC

    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
function calc_BS!(plank, proc, k_mtb, b_k_mtb, R_NC, R_PC)
    @inbounds proc.BS  .= min.(plank.Cq, plank.Nq ./ R_NC, plank.Pq ./ R_PC) .*
        k_mtb .* plank.Sz .^ b_k_mtb .* plank.ac
    @inbounds proc.exu .= max.(0.0, plank.Cq .- min.(plank.Cq, plank.Nq ./ R_NC, plank.Pq ./ R_PC)) .*
        k_mtb .* plank.Sz .^ b_k_mtb .* plank.ac

    return nothing
end

##### update C, N, P quotas, biomass, Chla, cell size
function update_biomass!(plank, proc, R_NC, R_PC, Cquota, Nsuper, ΔT)
    @inbounds plank.Bm  .= plank.Bm  .+ ΔT .*  proc.BS
    @inbounds plank.Cq  .= plank.Cq  .- ΔT .* (proc.BS .+ proc.exu)
    @inbounds plank.Nq  .= plank.Nq  .- ΔT .*  proc.BS .* R_NC
    @inbounds plank.Pq  .= plank.Pq  .- ΔT .*  proc.BS .* R_PC
    @inbounds plank.chl .= plank.chl .+ ΔT .*  proc.BS .* R_NC .* proc.ρchl
    @inbounds plank.age .= plank.age .+ ΔT ./ 3600.0 .* plank.ac
    @inbounds plank.Sz  .= (plank.Bm .+ plank.Cq) ./ (Cquota .* Nsuper)

    return nothing
end

##### calculate probability of cell division
##### sizer
function calc_dvid_size!(plank, proc, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    @inbounds proc.dvid .= dvid_P .* (tanh.(dvid_stp .* (plank.Sz .- dvid_reg)) .+ 1.0)
    @inbounds proc.dvid .= proc.dvid .* isless.(2*Cquota*Nsuper, plank.Bm) .* plank.ac

    return nothing
end
##### adder
function calc_dvid_add!(plank, proc, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    @inbounds proc.dvid .= dvid_P .* (tanh.(dvid_stp .* (plank.Sz .- plank.iS .- dvid_reg)) .+ 1.0)
    @inbounds proc.dvid .= proc.dvid .* isless.(2*Cquota*Nsuper, plank.Bm) .* plank.ac

    return nothing
end
##### age
function calc_dvid_age!(plank, proc, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper)
    @inbounds proc.dvid .= dvid_P .* (tanh.(dvid_stp .* (plank.age .- dvid_reg)) .+ 1.0)
    @inbounds proc.dvid .= proc.dvid .* isless.(2*Cquota*Nsuper, plank.Bm) .* plank.ac

    return nothing
end
##### timer
function calc_dvid_time!(plank, proc, dvid_stp, dvid_P, dvid_reg, Cquota, Nsuper, t)
    @inbounds proc.dvid .= dvid_P .* (tanh(dvid_stp * (t % 86400 ÷ 3600 - dvid_reg)) + 1.0)
    @inbounds proc.dvid .= proc.dvid .* isless.(2*Cquota*Nsuper, plank.Bm) .* plank.ac

    return nothing
end
##### timer & sizer
function calc_dvid_ts!(plank, proc, dvid_stp, dvid_stp2, dvid_P, dvid_reg, dvid_reg2, Cquota, Nsuper, t)
    @inbounds proc.dvid .= dvid_P .* (tanh(dvid_stp * (t % 86400 ÷ 3600 - dvid_reg)) + 1.0) .*
                                     (tanh.(dvid_stp2 .* (plank.Sz .- dvid_reg2)) .+ 1.0)
    @inbounds proc.dvid .= proc.dvid .* isless.(2*Cquota*Nsuper, plank.Bm) .* plank.ac

    return nothing
end

##### calculate the probability of grazing
##### quadratic grazing
function calc_graz_quadratic!(nuts, proc, grz_P)
    @inbounds proc.grz .= nuts.pop ./ grz_P

    return nothing
end
##### linear grazing decrease with depth
function calc_graz_linear!(plank, proc, grz_P, grz_stp)
    @inbounds proc.grz .= 1.0 ./ grz_P .* max.(0.15, 1.0 .- abs.(plank.x) ./ grz_stp)

    return nothing
end

##### calculate the probability of mortality
function calc_mort!(plank, proc, mort_reg, mort_P)
    @inbounds proc.mort .= mort_P .* (tanh.(6.0 .* (mort_reg .- plank.Sz)) .+ 1.0)

    return nothing
end

##### generate random numbers for grazing, mortality and division
function gen_rand_plk!(rnd, arch)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    return nothing
end

##### generate the random results from probabilities of grazing, mortality and cell division
function get_rands!(proc, rnd)
    @inbounds proc.grz  .= isless.(rnd.x, proc.grz)
    @inbounds proc.mort .= isless.(rnd.y, proc.mort)
    @inbounds proc.dvid .= isless.(rnd.z, proc.dvid)

    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    i = @index(Global, Linear)
    gi = @index(Group)
    @inbounds ctsdic[x[i], y[i], z[i], gi] += (proc.resp[i] - proc.PS[i]) * ΔT .* ac[i]
    @inbounds ctsdoc[x[i], y[i], z[i], gi] += (proc.exu[i] - proc.VDOC[i]) * ΔT .* ac[i]
    @inbounds ctsnh4[x[i], y[i], z[i], gi] -= proc.VNH4[i] * ΔT .* ac[i]
    @inbounds ctsno3[x[i], y[i], z[i], gi] -= proc.VNO3[i] * ΔT .* ac[i]
    @inbounds ctspo4[x[i], y[i], z[i], gi] -= proc.VPO4[i] * ΔT .* ac[i]
    # pay attention to ctspop[0,0,0,gi] for non-active individuals
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT, arch::Architecture)
    kernel! = calc_consume_kernel!(device(arch), 1, (size(ac,1),))
    event = kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                                   x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    gi = @index(Group)
    @inbounds ctsdoc[x[i], y[i], z[i], gi] += (plank.Bm[i] + plank.Cq[i]) * lossFracC .* plank.ac[i]
    @inbounds ctspoc[x[i], y[i], z[i], gi] += (plank.Bm[i] + plank.Cq[i]) * (1.0 - lossFracC) .* plank.ac[i]
    @inbounds ctsdon[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_NC + plank.Nq[i]) * lossFracN .* plank.ac[i]
    @inbounds ctspon[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_NC + plank.Nq[i]) * (1.0 - lossFracN) .* plank.ac[i]
    @inbounds ctsdop[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_PC + plank.Pq[i]) * lossFracP .* plank.ac[i]
    @inbounds ctspop[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_PC + plank.Pq[i]) * (1.0 - lossFracP) .* plank.ac[i]
    # pay attention to ctspop[0,0,0,gi] for non-active individuals
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                    x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC, arch::Architecture)
    kernel! = calc_loss_kernel!(device(arch), 1, (size(plank,1),))
    event = kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                    x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end
