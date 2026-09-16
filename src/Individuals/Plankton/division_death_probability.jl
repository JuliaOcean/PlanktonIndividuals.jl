abstract type division_type end
struct Sizer <: division_type end
struct Timer <: division_type end
struct Sizer_Timer <: division_type end

##### calculate probability of cell division
@inline calc_division(::Sizer, Sz, t, reg, reg2, P) = P * shape_func_inc(Sz, reg, 2.0f-3; pow = 2.0f0)

@inline calc_division(::Timer, Sz, t, reg, reg2, P) = P * shape_func_inc(t%86400/3600, reg, 1.0f-5)

@inline calc_division(::Sizer_Timer, Sz, t, reg, reg2, P) =
                        P * shape_func_inc(Sz, reg, 1.0f-3) * shape_func_inc(t%86400/3600, reg2, 1.0f-5)

@inline function divide_type(dvid_type)
    if dvid_type == 1
        return Sizer()
    elseif dvid_type == 2
        return Timer()
    elseif dvid_type == 3
        return Sizer_Timer()
    else
        throw(ArgumentError("Wrong cell division type, must be in 1 to 3"))
    end
end

@kernel function calc_dvid_kernel!(plank, dvid_type, p, t)
    i = @index(Global)
    @inbounds plank.dvid[i] = calc_division(dvid_type, plank.Sz[i], t,
                                            p.dvid_reg, p.dvid_reg2, p.dvid_P) * plank.ac[i]
end
function calc_dvid!(plank, dvid_type, p, t, arch)
    kernel! = calc_dvid_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, dvid_type, p, t)
    return nothing
end

@kernel function calc_MM_dvid_kernel!(plank, p)
    i = @index(Global)
    @inbounds C_struct = plank.PRO_RB[i] + plank.PRO_MC[i] + plank.PRO_MN[i] + 
                         plank.PRO_TN[i] + plank.PRO_TP[i] + plank.PRO_TFe[i] + 
                         plank.PRO_RS[i] + plank.PRO_PS[i] + plank.PRO_OT[i] +
                         plank.DNA[i] + plank.RNA[i] + plank.Chl[i] / 893.49f0 * 55.0f0

    @inbounds Qc = C_struct / (p.Cquota * p.Nsuper)
    @inbounds dvid_C =  shape_func_inc(Qc, p.dvid_reg, 1.0f-3; pow = 2.0f0) 
    @inbounds dvid_DNA = isless(2.0f0, plank.DNA[i] / (p.C_DNA * p.Nsuper)) 
    @inbounds plank.dvid[i] = p.dvid_P * dvid_C * dvid_DNA * plank.ac[i]
end
function calc_MM_dvid!(plank, p, arch)
    kernel! = calc_MM_dvid_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of grazing
##### quadratic grazing
@kernel function calc_graz_quadratic_kernel!(plank, trs, P)
    i = @index(Global)
    @inbounds plank.graz[i] = trs.pop[i] * P * plank.ac[i]
end
function calc_graz_quadratic!(plank, trs, P, arch)
    kernel! = calc_graz_quadratic_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, trs, P)
    return nothing
end

##### calculate the probability of mortality
@kernel function calc_mort_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.mort[i] = p.mort_P * shape_func_dec(plank.Sz[i], p.mort_reg, 1.0f-5) * plank.ac[i]
end
function calc_mort!(plank, p, arch)
    kernel! = calc_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of mortality caused by thermal exposure
##### when damaged biomass is greater than 99% of total biomass
@kernel function calc_thermal_mort_kernel!(plank, p)
    i = @index(Global)
    #@inbounds plank.mort[i] = p.mort_P * (1.0 - isless(plank.Bd[i], plank.Bm[i]*0.99)) * plank.ac[i]
    @inbounds plank.mort[i] = p.mort_P * shape_func_inc(plank.Bd[i]/max(1.0f-30, plank.Bm[i]), 0.99f0, 1.0f-6) * plank.ac[i]
end
function calc_thermal_mort!(plank, p, arch)
    kernel! = calc_thermal_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of mortality for MMM
@kernel function calc_MM_mort_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.mort[i] = p.mort_P * shape_func_inc(plank.age[i], p.mort_reg, 1f-4) * plank.ac[i]
end
function calc_MM_mort!(plank, p, arch)
    kernel! = calc_MM_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### generate the random results from probabilities of grazing, mortality and cell division
function get_probability!(plank, rnd, ΔT, arch)
    ##### generate random numbers (0,1) 
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    ##### compare the random number with the given probability
    ##### return 1 if random number is smaller
    @inbounds plank.graz .= isless.(rnd.x, plank.graz .* ΔT) .* plank.ac
    @inbounds plank.mort .= isless.(rnd.y, plank.mort .* ΔT) .* plank.ac
    @inbounds plank.dvid .= isless.(rnd.z, plank.dvid .* ΔT) .* plank.ac
    return nothing
end
