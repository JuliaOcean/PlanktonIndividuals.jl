##### temperature function
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

@kernel function calc_inorganic_uptake_kernel!(plank, proc, nuts, p)
    i = @index(Global)
    @inbounds proc.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.chl[i], plank.Bm[i], plank.Sz[i], p) * plank.ac[i]
end
function calc_inorganic_uptake!(plank, proc, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, nuts, p)
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

##### update C quotas
@kernel function update_quotas_kernel!(plank, proc, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += (proc.PS[i] - proc.resp[i]) * ΔT
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_quotas!(plank, proc, ΔT, arch)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, proc, ΔT)
    wait(device(arch), event)
    return nothing
end

##### update cell size
@kernel function update_cellsize_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.Sz[i]  = plank.Bm[i] / (p.Cquota * p.Nsuper)
    @inbounds plank.chl[i] = plank.Bm[i] * p.Chl2C
end
function update_cellsize!(plank, p, arch)
    kernel! = update_cellsize_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, p)
    wait(device(arch), event)
    return nothing
end