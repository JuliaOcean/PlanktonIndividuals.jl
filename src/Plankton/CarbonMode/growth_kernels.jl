##### temperature function
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp+273.15 - p.T⁺))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.T⁺-2)))
    return k/OGT_rate
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, temp, Chl, Bm, Sz, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * Sz^p.PC_b * tempFunc(temp, p)
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-10, PCm) * Chl / max(1.0e-10, Bm))) * Bm
    return PS
end

@kernel function calc_inorganic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.Chl[i], plank.Bm[i], plank.Sz[i], p) * plank.ac[i]

end

function calc_inorganic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, nuts, p)
    wait(device(arch), event)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.resp[i] = p.respir_a * plank.Sz[i]^p.respir_b * plank.Bm[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, T, p)
    wait(device(arch), event)
    return nothing
end

##### update C quotas
@kernel function update_quotas_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += (plank.PS[i] - plank.resp[i]) * ΔT
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_quotas!(plank, ΔT, arch)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, ΔT)
    wait(device(arch), event)
    return nothing
end

##### update cell size
@kernel function update_cellsize_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.Sz[i]  = plank.Bm[i] / (p.Cquota * p.Nsuper)
    @inbounds plank.Chl[i] = plank.Bm[i] * p.Chl2C
end
function update_cellsize!(plank, p, arch)
    kernel! = update_cellsize_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, p)
    wait(device(arch), event)
    return nothing
end

##### track temperature history
@kernel function calc_thermal_history_kernel!(plank, nuts, p, ΔT)
    i = @index(Global)
    @inbounds plank.Th[i] = plank.Th[i] * 0.9 + (nuts.T[i] + 273.15 - p.T⁺) * ΔT / 3600
end
function calc_thermal_history!(plank, nuts, p, ΔT, arch)
    kernel! = calc_thermal_history_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, nuts, p, ΔT)
    wait(device(arch), event)
    return nothing
end