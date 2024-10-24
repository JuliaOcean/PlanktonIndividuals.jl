##### temperature function for photosynthesis
@inline function tempFunc_PS(temp, p)
    x = temp - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 1.0f-1)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0))) * regT
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### temperature function for metabolic rates
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145f0*(temp+273.15f0)))
    k = max(0.0f0, k)
    OGT_rate = exp(-p.Ea/(8.3145f0*(p.Topt+273.15f0)))
    return min(1.0f0, k/OGT_rate)
end

##### allocation of functional biomass to repair
##### only functional when damaged biomass is greater than 0.0
##### use Bd/Bm as the allocation, more damaged biomass means more allocation to repair
@inline function gamma_alloc(Bm, Bd, p)
    γ = Bd / max(1.0f-30, Bm) * isless(0.0f0, Bd)
    γ = min(1.0f0, γ) * isequal(1.0f0, p.thermal)
    return γ
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_photosynthesis(par, temp, Chl, Bm, Bd, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * (tempFunc(temp, p) * p.thermal + tempFunc_PS(temp, p) * (1.0f0 - p.thermal))
    light_limit = 1.0f0 - exp(-αI / max(1.0f-30, PCm) * Chl / max(1.0f-30, Bm))
    PS  = PCm * max(0.0f0, Bm - Bd) * (light_limit * (1.0f0 - p.is_bact) + p.is_bact)
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
@kernel function calc_repair_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.RP[i] = plank.PS[i] * gamma_alloc(plank.Bm[i], plank.Bd[i], p)
end
function calc_repair!(plank, p, arch)
    kernel! = calc_repair_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate biosynthesis rate (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.BS[i] = plank.PS[i] * (1.0f0 - gamma_alloc(plank.Bm[i], plank.Bd[i], p))
end
function calc_BS!(plank, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate thermal damage rate (mmolC/individual/second)
##### keep thermal damage between 0.0 and healthy biomass (Bm-Bd)
@kernel function calc_thermal_damage_kernel!(plank, T, p, ΔT)
    i = @index(Global)
    @inbounds plank.TD[i] = (T[i] - p.Topt) * p.f_T2B * (plank.Bm[i] - plank.Bd[i]) *
                            isless(p.Topt, T[i]) *
                            isless(plank.Bd[i], plank.Bm[i]) *
                            isequal(1.0f0, p.thermal)
    @inbounds plank.TD[i] = max(0.0f0, min(plank.TD[i], (plank.Bm[i] - plank.Bd[i]) / ΔT))
end
function calc_thermal_damage!(plank, T, p, ΔT, arch)
    kernel! = calc_thermal_damage_kernel!(device(arch), 256, (size(plank.ac, 1)))
    kernel!(plank, T, p, ΔT)
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
    @inbounds plank.Bm[i]  += (plank.BS[i] - plank.RS[i]) * ΔT
    @inbounds plank.Bd[i]  += (plank.TD[i] - plank.RP[i]) * ΔT
    @inbounds plank.age[i] += ΔT / 3600.0f0 * plank.ac[i]
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
