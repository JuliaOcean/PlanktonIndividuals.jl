##### temperature function
# @inline function tempFunc(temp, p)
#     TFunc = max(1.0e-10, exp(p.TAe * (1.0/(temp+273.15) - (1.0/p.TRef)))) * p.TCoeff
#     return TFunc
# end
@inline function tempFunc(temp, p)
    TFunc = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp+273.15 - p.T⁺))
    return max(0.0, TFunc)
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, temp, chl, Bm, Sz, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * Sz^p.PC_b * tempFunc(temp, p)
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-10, PCm) * chl / max(1.0e-10, Bm))) * Bm
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NH4_uptake(NH4, temp, Cq, Nq, Bm, Sz, p)
    regQ = max(0.0, min(1.0, (p.Nqmax - (Nq + Bm * p.R_NC) / max(1.0e-10, Bm + Cq)) / (p.Nqmax - p.Nqmin)))
    VN = p.VNH4max * Sz^p.VN_b * regQ * NH4/max(1.0e-10, NH4+p.KsatNH4) * tempFunc(temp, p) * Bm
    return VN
end
@inline function calc_NO3_uptake(NO3, temp, Cq, Nq, Bm, Sz, p)
    regQ = max(0.0, min(1.0, (p.Nqmax - (Nq + Bm * p.R_NC) / max(1.0e-10, Bm + Cq)) / (p.Nqmax - p.Nqmin)))
    VN = p.VNO3max * Sz^p.VN_b * regQ * NO3/max(1.0e-10, NO3+p.KsatNO3) * tempFunc(temp, p) * Bm
    return VN
end
@inline function calc_PO4_uptake(PO4, temp, Cq, Pq, Bm, Sz, p)
    regQ = max(0.0, min(1.0, (p.Pqmax - (Pq + Bm * p.R_PC) / max(1.0e-10, Bm + Cq)) / (p.Pqmax - p.Pqmin)))
    VN = p.VPO4max * Sz^p.VP_b * regQ * PO4/max(1.0e-10, PO4+p.KsatPO4) * tempFunc(temp, p) * Bm
    return VN
end
@inline function calc_DOC_uptake(DOC, temp, Cq, Bm, Sz, p)
    regQ = max(0.0, min(1.0, (p.Cqmax - Cq / max(1.0e-10, Bm + Cq)) / (p.Cqmax - p.Cqmin)))
    VN = p.VDOCmax * Sz^p.VDOC_b * regQ * DOC/max(1.0e-10, DOC+p.KsatDOC) * tempFunc(temp, p) * Bm
    return VN
end

@kernel function calc_inorganic_uptake_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.chl[i], plank.Bm[i], plank.Sz[i], p) * plank.ac[i]

    @inbounds proc.VNH4[i] = calc_NH4_uptake(nuts.NH4[i], nuts.T[i], plank.Cq[i], plank.Nq[i], plank.Bm[i],
                                             plank.Sz[i], p) * plank.ac[i]

    @inbounds proc.VNO3[i] = calc_NO3_uptake(nuts.NO3[i], nuts.T[i], plank.Cq[i], plank.Nq[i], plank.Bm[i],
                                             plank.Sz[i], p) * plank.ac[i]

    @inbounds proc.VPO4[i] = calc_PO4_uptake(nuts.PO4[i], nuts.T[i], plank.Cq[i], plank.Pq[i], plank.Bm[i],
                                             plank.Sz[i], p) * plank.ac[i]
end
function calc_inorganic_uptake!(plank, proc, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, p)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas for the first time of each time step
@kernel function update_quotas_1_kernel!(plank, proc, ΔT)
    i = @index(Global)
    @inbounds plank.Cq[i] += proc.PS[i]   * ΔT
    @inbounds plank.Pq[i] += proc.VPO4[i] * ΔT
    @inbounds plank.Nq[i] +=(proc.VNH4[i] + proc.VNO3[i]) * ΔT
end
function update_quotas_1!(plank, proc, ΔT, arch)
    kernel! = update_quotas_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, ΔT)
    wait(device(arch), event)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 5% of total C acquisition.
@kernel function calc_organic_uptake_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.VDOC[i] = calc_DOC_uptake(nuts.DOC[i], nuts.T[i], plank.Cq[i], plank.Bm[i], plank.Sz[i], p) * plank.ac[i]

    @inbounds proc.VDOC[i] = proc.VDOC[i] * isless(0.05, proc.PS[i]/(proc.VDOC[i]+proc.PS[i]))
end
function calc_organic_uptake!(plank, proc, nuts, p, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, p)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(plank, proc, par, p)
    i = @index(Global)
    @inbounds proc.ρchl[i] = proc.PS[i]/max(1.0e-10, plank.Bm[i]) * p.Chl2N / 
                             max(1.0e-10, par[i] * p.α * p.Φ * plank.chl[i]/plank.Bm[i]) *
                             isless(1.0e-1, par[i]) * plank.ac[i]
end
function calc_ρchl!(plank, proc, par, p, arch)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, par, p)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, proc, T, p)
    i = @index(Global)
    @inbounds proc.resp[i] = p.respir_a * plank.Sz[i]^p.respir_b * plank.Bm[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, proc, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, T, p)
    wait(device(arch), event)
    return nothing
end


##### update C, N, P quotas for the second time of each time step
@kernel function update_quotas_2_kernel!(plank, proc, ΔT, p)
    i = @index(Global)
    @inbounds plank.Cq[i] = plank.Cq[i] + (proc.VDOC[i] - proc.resp[i]) * ΔT
    @inbounds plank.Cq[i] = plank.Cq[i] + max(0.0, (0.0 - plank.Cq[i]))
    @inbounds plank.Bm[i] = plank.Bm[i] - max(0.0, (0.0 - plank.Cq[i]))
    @inbounds plank.Nq[i] = plank.Nq[i] + max(0.0, (0.0 - plank.Cq[i])) * p.R_NC
    @inbounds plank.Pq[i] = plank.Pq[i] + max(0.0, (0.0 - plank.Cq[i])) * p.R_PC
end
function update_quotas_2!(plank, proc, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, ΔT, p)
    wait(device(arch), event)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, proc, p)
    i = @index(Global)
    @inbounds proc.BS[i] = min(plank.Cq[i], plank.Nq[i]/p.R_NC, plank.Pq[i]/p.R_PC) * 
                           p.k_mtb * plank.Sz[i]^p.k_mtb_b * plank.ac[i]
    @inbounds proc.exu[i]= max(0.0, plank.Cq[i] - min(plank.Cq[i], plank.Nq[i]/p.R_NC, plank.Pq[i]/p.R_PC)) *
                           p.k_mtb * plank.Sz[i]^p.k_mtb_b * plank.ac[i]
end
function calc_BS!(plank, proc, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, p)
    wait(device(arch), event)
    return nothing
end

##### update C, N, P quotas, biomass, Chla
@kernel function update_biomass_kernel!(plank, proc, p, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += ΔT * proc.BS[i]
    @inbounds plank.Cq[i]  -= ΔT *(proc.BS[i] + proc.exu[i])
    @inbounds plank.Nq[i]  -= ΔT * proc.BS[i] * p.R_NC
    @inbounds plank.Pq[i]  -= ΔT * proc.BS[i] * p.R_PC
    @inbounds plank.chl[i] += ΔT * proc.BS[i] * p.R_NC * proc.ρchl[i]
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_biomass!(plank, proc, p, ΔT, arch)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, p, ΔT)
    wait(device(arch), event)
    return nothing
end

##### update cell size
@kernel function update_cellsize_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.Sz[i] = (plank.Bm[i] + plank.Cq[i]) / (p.Cquota * p.Nsuper)
end
function update_cellsize!(plank, p, arch)
    kernel! = update_cellsize_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, p)
    wait(device(arch), event)
    return nothing
end