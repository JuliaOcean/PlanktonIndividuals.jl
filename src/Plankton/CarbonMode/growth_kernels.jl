##### temperature function for photosynthesis
@inline function tempFunc_PS(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp - p.Tmax))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.Topt+273.15)))
    return min(1.0, k/OGT_rate)
end

##### temperature function for metabolic rates
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.Topt+273.15)))
    return min(1.0, k/OGT_rate)
end

##### allocation of functional biomass between biosynthesis and repair
@inline function gamma_alloc(Bm, Bd)
    # γ = max(0.0, temp - p.Topt) / (p.Tmax - p.Topt)
    γ = Bd / max(1.0e-30, Bm)
    γ = min(1.0, γ)
    return γ
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_photosynthesis(par, temp, Chl, Bm, Bd, p)
    αI  = par * p.α * p.Φ
    if p.thermal == 1.0
        PCm = p.PCmax * tempFunc(temp, p)
    else
        PCm = p.PCmax * tempFunc_PS(temp, p)
    end
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-30, PCm) * Chl / max(1.0e-30, Bm))) * max(0.0, Bm - Bd)
    return PS
end

@kernel function calc_PS_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_photosynthesis(nuts.par[i], nuts.T[i], plank.Chl[i], plank.Bm[i], plank.Bd[i], p) * plank.ac[i]
end

function calc_PS!(plank, nuts, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### calculate repair rate (mmolC/individual/second)
@kernel function calc_repair_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.RP[i] = plank.PS[i] * gamma_alloc(plank.Bm[i], plank.Bd[i])
end
function calc_repair!(plank, T, p, arch)
    kernel! = calc_repair_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate biosynthesis rate (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.BS[i] = plank.PS[i] * (1.0 - gamma_alloc(plank.Bm[i], plank.Bd[i]))
end
function calc_BS!(plank, T, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate thermal damage rate (mmolC/individual/second)
@kernel function calc_thermal_damage_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.TD[i] = (T[i] - p.Topt) * p.f_T2B * isless(p.Topt, T[i]) * isless(0.0, p.thermal)
end
function calc_thermal_damage!(plank, T, p, arch)
    kernel! = calc_thermal_damage_kernel!(device(arch), 256, (size(plank.ac, 1)))
    kernel!(plank, T, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.RS[i] = p.respir * plank.Bm[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end

##### update C quotas
@kernel function update_quotas_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += (plank.BS[i] - plank.RS[i] - plank.TD[i] + plank.RP[i]) * ΔT
    @inbounds plank.Bd[i]  += (plank.TD[i] - plank.RP[i]) * ΔT
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_quotas!(plank, ΔT, arch)
    kernel! = update_quotas_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
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
    kernel!(plank, p)
    return nothing
end
