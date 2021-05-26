##### calculate Chla and individual counts based on the status of plankton individuals
@inline function gpu_acc_counts_kernel!(ctschl, ctspop, chl, ac, x, y, z)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds @atomic ctschl[x[i], y[i], z[i]] += chl[i] * ac[i]
        @inbounds @atomic ctspop[x[i], y[i], z[i]] += ac[i]
    end
    return nothing
end
function acc_counts!(ctschl, ctspop, chl, ac, x, y, z, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_acc_counts_kernel!(ctschl, ctspop, chl, ac, x, y, z)
    return nothing 
end
function acc_counts!(ctschl, ctspop, chl, ac, x, y, z, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctschl[x[i], y[i], z[i]] += chl[i] * ac[i]
        @inbounds ctspop[x[i], y[i], z[i]] += ac[i]
    end
    return nothing
end

##### deal with nutrients uptake
function gpu_calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds @atomic ctsdic[x[i], y[i], z[i]] += (proc.resp[i] - proc.PS[i]) * ΔT * ac[i]
        @inbounds @atomic ctsdoc[x[i], y[i], z[i]] += (proc.exu[i] - proc.VDOC[i]) * ΔT * ac[i]
        @inbounds @atomic ctsnh4[x[i], y[i], z[i]] -= proc.VNH4[i] * ΔT * ac[i]
        @inbounds @atomic ctsno3[x[i], y[i], z[i]] -= proc.VNO3[i] * ΔT * ac[i]
        @inbounds @atomic ctspo4[x[i], y[i], z[i]] -= proc.VPO4[i] * ΔT * ac[i]
    end
    return nothing
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, 
                                                        ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    return nothing 
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdic[x[i], y[i], z[i]] += (proc.resp[i] - proc.PS[i]) * ΔT * ac[i]
        @inbounds ctsdoc[x[i], y[i], z[i]] += (proc.exu[i] - proc.VDOC[i]) * ΔT * ac[i]
        @inbounds ctsnh4[x[i], y[i], z[i]] -= proc.VNH4[i] * ΔT * ac[i]
        @inbounds ctsno3[x[i], y[i], z[i]] -= proc.VNO3[i] * ΔT * ac[i]
        @inbounds ctspo4[x[i], y[i], z[i]] -= proc.VPO4[i] * ΔT * ac[i]
    end
    return nothing
end

##### deal with grazed or dead individuals
function gpu_calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                                   loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds @atomic ctsdoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * lossFracC * ac[i] * loss[i]
        @inbounds @atomic ctsdon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * lossFracN * ac[i] * loss[i]
        @inbounds @atomic ctsdop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * lossFracP * ac[i] * loss[i]
        @inbounds @atomic ctspoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * (1.0-lossFracC) * ac[i] * loss[i]
        @inbounds @atomic ctspon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * (1.0-lossFracN) * ac[i] * loss[i]
        @inbounds @atomic ctspop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * (1.0-lossFracP) * ac[i] * loss[i]
    end
    return nothing
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                    loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, 
                                                        ctsdop, ctspop, plank, ac, x, y, z, loss, 
                                                        lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    return nothing 
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                           loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * lossFracC * ac[i] * loss[i]
        @inbounds ctsdon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * lossFracN * ac[i] * loss[i]
        @inbounds ctsdop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * lossFracP * ac[i] * loss[i]
        @inbounds ctspoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * (1.0-lossFracC) * ac[i] * loss[i]
        @inbounds ctspon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * (1.0-lossFracN) * ac[i] * loss[i]
        @inbounds ctspop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * (1.0-lossFracP) * ac[i] * loss[i]
    end
    return nothing
end

using KernelAbstractions.Extras.LoopInfo: @unroll
##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, chl, PARF, g::AbstractGrid, kc, kw)
    i, j = @index(Global, NTuple)
    @unroll for k in 1:g.Nz
        ii = i + g.Hx
        jj = j + g.Hy
        kk = k + g.Hz
        atten = (chl[ii,jj,kk]/volume(ii, jj, kk, g) * kc + kw) * g.Δz
        par[ii,jj,kk] = PARF[i,j] * (1.0 - exp(-atten)) / atten
        PARF[i,j] = PARF[i,j] * exp(-atten)
    end
end
function calc_par!(par, arch::Architecture, chl, PARF, g::AbstractGrid, kc, kw)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    event = kernel!(par, chl, PARF, g, kc, kw)
    wait(device(arch), event)
    return nothing
end