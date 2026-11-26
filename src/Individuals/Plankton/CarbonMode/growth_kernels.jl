##### temperature function for photosynthesis
@inline function tempFunc_PS(T, p)
    x = T - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 1.0f-1)
    Ea = p.Ea_ref * (p.Topt_ref / max(1.0f-30, p.Topt))
    Δtemp = 1.0f0/(T+273.15f0) - 1.0f0/(p.Topt+273.15f0)
    k = exp(-Ea/8.3145f0*Δtemp) * regT
    return min(1.0f0, k)
end

##### temperature function for metabolic rates
@inline function tempFunc(T, p)
    Ea = p.Ea_ref * (p.Topt_ref / p.Topt)
    Δtemp = 1.0f0/(T+273.15f0) - 1.0f0/(p.Topt+273.15f0)
    k = exp(-Ea/8.3145f0*Δtemp)
    return min(1.0f0, k)
end

@kernel function calc_nut_limit_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.fNO3[i] = trs.NO3[i] / max(1.0f-30, p.KsatNO3 + trs.NO3[i])
    @inbounds plank.fNH4[i] = trs.NH4[i] / max(1.0f-30, p.KsatNH4 + trs.NH4[i])
    @inbounds plank.fPO4[i] = trs.PO4[i] / max(1.0f-30, p.KsatPO4 + trs.PO4[i])
end
function calc_nut_limit!(plank, trs, p, arch)
    kernel! = calc_nut_limit_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### allocation of functional biomass to growth and repair
##### only functional when damaged biomass is greater than 0.0
##### use Bd/Bm as the allocation factor
##### more damaged biomass means more allocation to repair
@inline function gamma_alloc(Bm, Bd, p)
    γ = Bd / max(1.0f-30, Bm)
    γ = min(1.0f0, γ) * isequal(1.0f0, p.thermal)
    return 1.0f0 - γ
end
@kernel function gamma_alloc_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.γ[i] = gamma_alloc(plank.Bm[i], plank.Bd[i], p)
    @inbounds plank.γ[i] = plank.γ[i] * min(plank.fNO3[i] + plank.fNH4[i], plank.fPO4[i])
end
function calc_gamma_alloc!(plank, p, arch)
    kernel! = gamma_alloc_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_photosynthesis(par, T, Chl, Bm, Bd, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * (tempFunc(T, p) * p.thermal + tempFunc_PS(T, p) * (1.0f0 - p.thermal))
    light_limit = 1.0f0 - exp(-αI / max(1.0f-30, PCm) * Chl / max(1.0f-30, Bm - Bd))
    PS  = PCm * max(0.0f0, Bm - Bd) * (light_limit * (1.0f0 - p.is_bact) + p.is_bact)
    return PS
end

@kernel function calc_PS_kernel!(plank, trs, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_photosynthesis(trs.par[i], trs.T[i], plank.Chl[i],
                            plank.Bm[i], plank.Bd[i], p) * plank.ac[i]
end

function calc_PS!(plank, trs, p, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, p)
    return nothing
end

##### calculate repair rate (mmolC/individual/second)
@kernel function calc_C_alloc_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.BS[i] = plank.PS[i] * plank.γ[i]
    @inbounds plank.RP[i] = plank.PS[i] * (1.0f0 - plank.γ[i])
    @inbounds plank.RP[i] = max(0.0f0, min(plank.RP[i], plank.Bd[i] / ΔT))
    @inbounds plank.RS[i] = plank.PS[i] - plank.RP[i] - plank.BS[i]
end
function calc_C_alloc!(plank, p, ΔT, arch)
    kernel! = calc_C_alloc_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
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

##### update C quotas
@kernel function update_quotas_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += plank.BS[i] * ΔT
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
