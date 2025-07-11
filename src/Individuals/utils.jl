##### find indices (halo points included)
@kernel function find_inds_kernel!(particle, g::AbstractGrid)
    i = @index(Global)
    @inbounds particle.xi[i] = unsafe_trunc(Int, get_xf_index(particle.x[i]) * particle.ac[i]) + g.Hx 
    @inbounds particle.yi[i] = unsafe_trunc(Int, get_yf_index(particle.y[i]) * particle.ac[i]) + g.Hy
    @inbounds particle.zi[i] = unsafe_trunc(Int, get_zf_index(particle.z[i]) * particle.ac[i]) + g.Hz
end
function find_inds!(particle, g::AbstractGrid, arch::Architecture)
    kernel! = find_inds_kernel!(device(arch), 256, (size(particle.ac,1)))
    kernel!(particle, g)
    return nothing
end

@kernel function find_NPT_kernel!(trs, x, y, z, ac, NH4, NO3, PO4, DOC, FeT, par, par₀, temp, pop)
    i = @index(Global)
    @inbounds trs.NH4[i] = max(0.0f0, NH4[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.NO3[i] = max(0.0f0, NO3[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.PO4[i] = max(0.0f0, PO4[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.DOC[i] = max(0.0f0, DOC[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.FeT[i] = max(0.0f0, FeT[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.par[i] = par[x[i], y[i], z[i]] * ac[i]
    @inbounds trs.dpar[i]= (par[x[i], y[i], z[i]] - par₀[x[i], y[i], z[i]]) * ac[i]
    @inbounds trs.T[i]   =temp[x[i], y[i], z[i]] * ac[i]
    @inbounds trs.pop[i] = pop[x[i], y[i], z[i]] * ac[i]
end
function find_NPT!(trs, x, y, z, ac, NH4, NO3, PO4, DOC, FeT, par, par₀, temp, pop, arch::Architecture)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(trs, x, y, z, ac, NH4, NO3, PO4, DOC, FeT, par, par₀, temp, pop)
    return nothing
end

##### calculate Chla based on the status of plankton individuals
@kernel function acc_chl_kernel!(ctsChl, Chl, ac, x, y, z)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsChl[x[i], y[i], z[i]] += Chl[i] * ac[i]
end
function acc_chl!(ctsChl, Chl, ac, x, y, z, arch)
    kernel! = acc_chl_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsChl, Chl, ac, x, y, z)
    return nothing 
end

##### calculate individual counts based on the status of plankton individuals
@kernel function acc_counts_kernel!(ctspop, ac, x, y, z)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctspop[x[i], y[i], z[i]] += ac[i]
end
function acc_counts!(ctspop, ac, x, y, z, arch)
    kernel! = acc_counts_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctspop, ac, x, y, z)
    return nothing 
end

##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, Chl, PARF, g::AbstractGrid, kc, kw, ki)
    i, j = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = ki + g.Hz
    atten = (Chl[ii,jj,kk]/volume(ii, jj, kk, g) * kc + kw) * ΔzF(ii, jj, kk, g)
    par[ii,jj,kk] = PARF[i,j] * (1.0f0 - exp(-atten)) / atten
    PARF[i,j] = PARF[i,j] * exp(-atten)
end
function calc_par!(par, arch::Architecture, Chl, PARF, g::AbstractGrid, kc, kw, ki)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    kernel!(par, Chl, PARF, g, kc, kw, ki)
    return nothing
end

##### mask individuals due to land shape
@kernel function mask_individuals_kernel!(particle, g::AbstractGrid)
    i = @index(Global)
    @inbounds xi = unsafe_trunc(Int, (particle.x[i]+1) * particle.ac[i]) + g.Hx 
    @inbounds yi = unsafe_trunc(Int, (particle.y[i]+1) * particle.ac[i]) + g.Hy
    @inbounds zi = unsafe_trunc(Int, (particle.z[i]+1) * particle.ac[i]) + g.Hz
    @inbounds particle.ac[i] = g.landmask[xi, yi, zi] * particle.ac[i]
end
function mask_individuals!(particle, g::AbstractGrid, N, arch)
    kernel! = mask_individuals_kernel!(device(arch), 256, (N,))
    kernel!(particle, g)
    return nothing
end

##### shape function - decrease from 1.0 to 0.0 while x increase from 0.0 to 1.0
##### sharp decrease near x/xmax = 1.0
@inline function shape_func_dec(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, 1.0f0 - x / xmax))
    reg = fx^pow / (k + fx^pow)
    return reg
end

##### shape function - increase from 0.0 to 1.0 while x increase from 0.0 to 1.0
##### sharp increase near x/xmax = 1.0
@inline function shape_func_inc(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, 1.0f0 - x / xmax))
    reg = fx^pow / (k + fx^pow)
    return 1.0f0 - reg
end

##### shape function - decrease from 1.0 to 0.0 while x increase from 0.0 to 1.0
##### sharp decrease near x/xmax = 0.0
@inline function shape_func_dec_alt(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, x / xmax))
    reg = fx^pow / (k + fx^pow)
    return 1.0f0 - reg
end

##### shape function - increase from 0.0 to 1.0 while x increase from 0.0 to 1.0
##### sharp increase near x/xmax = 0.0
@inline function shape_func_inc_alt(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, x / xmax))
    reg = fx^pow / (k + fx^pow)
    return reg
end