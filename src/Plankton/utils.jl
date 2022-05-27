##### structs for individuals
mutable struct plankton
    data::AbstractArray
    p::NamedTuple
end

struct individuals
    phytos::NamedTuple
    # zoos::NamedTuple
end

##### find indices (halo points included)
@kernel function find_inds_kernel!(plank, g::AbstractGrid)
    i = @index(Global)
    @inbounds plank.xi[i] = unsafe_trunc(Int, get_xf_index(plank.x[i]) * plank.ac[i]) + g.Hx 
    @inbounds plank.yi[i] = unsafe_trunc(Int, get_yf_index(plank.y[i]) * plank.ac[i]) + g.Hy
    @inbounds plank.zi[i] = unsafe_trunc(Int, get_zf_index(plank.z[i]) * plank.ac[i]) + g.Hz
end
function find_inds!(plank, g::AbstractGrid, arch::Architecture)
    kernel! = find_inds_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, g)
    wait(device(arch), event)
    return nothing
end

@kernel function find_NPT_kernel!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop)
    i = @index(Global)
    @inbounds nuts.NH4[i] = max(1.0e-30, NH4[x[i], y[i], z[i]]) * ac[i]
    @inbounds nuts.NO3[i] = max(1.0e-30, NO3[x[i], y[i], z[i]]) * ac[i]
    @inbounds nuts.PO4[i] = max(1.0e-30, PO4[x[i], y[i], z[i]]) * ac[i]
    @inbounds nuts.DOC[i] = max(1.0e-30, DOC[x[i], y[i], z[i]]) * ac[i]
    @inbounds nuts.par[i] = par[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.T[i]   =temp[x[i], y[i], z[i]] * ac[i]
    @inbounds nuts.pop[i] = pop[x[i], y[i], z[i]] * ac[i]
end
function find_NPT!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop, arch::Architecture)
    kernel! = find_NPT_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(nuts, x, y, z, ac, NH4, NO3, PO4, DOC, par, temp, pop)
    wait(device(arch), event)
    return nothing
end

##### calculate Chla and individual counts based on the status of plankton individuals
@inline function gpu_acc_counts_kernel!(ctsChl, ctspop, Chl, ac, x, y, z)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic ctsChl[x[i], y[i], z[i]] += Chl[i] * ac[i]
        @inbounds CUDA.@atomic ctspop[x[i], y[i], z[i]] += ac[i]
    end
    return nothing
end
function acc_counts!(ctsChl, ctspop, Chl, ac, x, y, z, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_acc_counts_kernel!(ctsChl, ctspop, Chl, ac, x, y, z)
    return nothing 
end
function acc_counts!(ctsChl, ctspop, Chl, ac, x, y, z, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsChl[x[i], y[i], z[i]] += Chl[i] * ac[i]
        @inbounds ctspop[x[i], y[i], z[i]] += ac[i]
    end
    return nothing
end

##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, Chl, PARF, g::AbstractGrid, kc, kw)
    i, j = @index(Global, NTuple)
    @unroll for k in 1:g.Nz
        ii = i + g.Hx
        jj = j + g.Hy
        kk = k + g.Hz
        atten = (Chl[ii,jj,kk]/volume(ii, jj, kk, g) * kc + kw) * ΔzF(ii, jj, kk, g)
        par[ii,jj,kk] = PARF[i,j] * (1.0 - exp(-atten)) / atten
        PARF[i,j] = PARF[i,j] * exp(-atten)
    end
end
function calc_par!(par, arch::Architecture, Chl, PARF, g::AbstractGrid, kc, kw)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    event = kernel!(par, Chl, PARF, g, kc, kw)
    wait(device(arch), event)
    return nothing
end

@kernel function mask_individuals_kernel!(plank, g::AbstractGrid)
    i = @index(Global)
    @inbounds xi = unsafe_trunc(Int, (plank.x[i]+1) * plank.ac[i]) + g.Hx 
    @inbounds yi = unsafe_trunc(Int, (plank.y[i]+1) * plank.ac[i]) + g.Hy
    @inbounds zi = unsafe_trunc(Int, (plank.z[i]+1) * plank.ac[i]) + g.Hz
    @inbounds plank.ac[i] = g.landmask[xi, yi, zi] * plank.ac[i]
end
function mask_individuals!(plank, g::AbstractGrid, N, arch)
    kernel! = mask_individuals_kernel!(device(arch), 256, (N,))
    event = kernel!(plank, g)
    wait(device(arch), event)
    return nothing
end