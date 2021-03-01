##### find indices (halo points included)
@kernel function find_inds_kernel!(plank, g::Grids)
    i = @index(Global)
    @inbounds plank.xi[i] = unsafe_trunc(Int, get_xf_index(plank.x[i], g) * plank.ac[i]) + g.Hx + 1
    @inbounds plank.yi[i] = unsafe_trunc(Int, get_yf_index(plank.y[i], g) * plank.ac[i]) + g.Hy + 1
    @inbounds plank.zi[i] = unsafe_trunc(Int, get_zf_index(plank.z[i], g) * plank.ac[i]) + g.Hz + 1
end
function find_inds!(plank, g::Grids, arch::Architecture)
    kernel! = find_inds_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, g)
    wait(device(arch), event)
    return nothing
end

@kernel function find_NPT_kernel!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop)
    i = @index(Global)
    @inbounds nuts.NH4[i] = NH4[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.NO3[i] = NO3[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.PO4[i] = PO4[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.DOC[i] = DOC[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.αI[i]  = par[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.Tem[i] =temp[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.pop[i] = pop[x[i], y[i], z[i]] * ac[i]
end
function find_NPT!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop, arch::Architecture)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop)
    wait(device(arch), event)
    return nothing
end

@kernel function update_NPT_kernel!(nuts, ac, p)
    i = @index(Global)
    @inbounds nuts.NH4[i] = max(1.0e-10, nuts.NH4[i]) * ac[i]
    @inbounds nuts.NO3[i] = max(1.0e-10, nuts.NO3[i]) * ac[i]
    @inbounds nuts.PO4[i] = max(1.0e-10, nuts.PO4[i]) * ac[i]
    @inbounds nuts.DOC[i] = max(1.0e-10, nuts.DOC[i]) * ac[i]
    @inbounds nuts.αI[i]  = nuts.αI[i] * p.α * p.Φ *ac[i]
    @inbounds nuts.Tem[i] = max(1.0e-10, exp(p.TempAe * (1.0/(nuts.Tem[i]+273.15) - (1.0/p.Tempref)))) * p.TempCoeff * ac[i]
end
function update_NPT!(nuts, ac, p, arch::Architecture)
    kernel! = update_NPT_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(nuts, ac, p)
    wait(device(arch), event)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(αI, temp, chl, Bm, Sz, PCmax, PC_b)
    PCm = PCmax * Sz^PC_b * temp
    PS  = PCm * (1.0 - exp(-αI * chl / max(1.0e-10, Bm * PCm))) * Bm
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_uptake(nut, temp, Cq, Nq, Bm, Sz, Nqmax, Nqmin, Vmax, VN_b, Ksat, R_NC)
    regQ = max(0.0, min(1.0, (Nqmax - (Nq + Bm * R_NC) / max(1.0e-10, Bm + Cq)) / (Nqmax - Nqmin)))
    VN = Vmax * Sz^VN_b * regQ * nut/max(1.0e-10, nut+Ksat) * temp * Bm
    return VN
end

@kernel function calc_inorganic_uptake_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.PS[i] = calc_PS(nuts.αI[i], nuts.Tem[i], plank.chl[i], plank.Bm[i], 
                                   plank.Sz[i], p.PCmax, p.PC_b) * plank.ac[i]

    @inbounds proc.VNH4[i] = calc_uptake(nuts.NH4[i], nuts.Tem[i], plank.Cq[i], plank.Nq[i], plank.Bm[i], plank.Sz[i],
                                         p.Nqmax, p.Nqmin, p.VNH4max, p.VN_b, p.KsatNH4, p.R_NC) * plank.ac[i]

    @inbounds proc.VNO3[i] = calc_uptake(nuts.NO3[i], nuts.Tem[i], plank.Cq[i], plank.Nq[i], plank.Bm[i], plank.Sz[i],
                                         p.Nqmax, p.Nqmin, p.VNO3max, p.VN_b, p.KsatNO3, p.R_NC) * plank.ac[i]

    @inbounds proc.VPO4[i] = calc_uptake(nuts.PO4[i], nuts.Tem[i], plank.Cq[i], plank.Pq[i], plank.Bm[i], plank.Sz[i],
                                         p.Pqmax, p.Pqmin, p.VPO4max, p.VN_b, p.KsatPO4, p.R_PC) * plank.ac[i]
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
@kernel function calc_organic_uptake_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.VDOC[i] = calc_uptake(nuts.DOC[i], nuts.Tem[i], plank.Cq[i], plank.Cq[i], plank.Bm[i], plank.Sz[i],
                                         p.Cqmax, p.Cqmin, p.VDOCmax, p.VDOC_b, p.KsatDOC, 0.0) * plank.ac[i]
end
function calc_organic_uptake!(plank, proc, nuts, p, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, p)
    wait(device(arch), event)
    return nothing
end

##### calculate ρchl
@kernel function calc_ρchl_kernel!(plank, proc, nuts, Chl2N)
    i = @index(Global)
    @inbounds proc.ρchl[i] = proc.PS[i]/max(1.0e-10, plank.Bm[i]) * Chl2N/max(1.0e-10, nuts.αI[i] * plank.chl[i]/plank.Bm[i]) *
                   isless(1.0e-8, nuts.αI[i]) * plank.ac[i]
end
function calc_ρchl!(plank, proc, nuts, Chl2N, arch)
    kernel! = calc_ρchl_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, Chl2N)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.resp[i] = p.respir_a * plank.Sz[i]^p.respir_b * plank.Bm[i] * nuts.Tem[i] * plank.ac[i]
end
function calc_respir!(plank, proc, nuts, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, p)
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