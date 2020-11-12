##### calculate Chla and individual counts based on the status of plankton individuals
@kernel function acc_counts_kernel!(ctschl, ctspop, chl, ac, x, y, z)
    i = @index(Global)
    gi = @index(Group)
    @inbounds ctschl[x[i], y[i], z[i], gi] += copy(chl[i]) * ac[i]
    @inbounds ctspop[x[i], y[i], z[i], gi] += 1.0 * ac[i]
    # pay attention to cumpop[0,0,0,gi] for non-active individuals
end
function acc_counts!(ctschl, ctspop, chl, ac, x, y, z, arch::Architecture)
    kernel! = acc_counts_kernel!(device(arch), 1, (size(ac,1),))
    event = kernel!(ctschl, ctspop, chl, ac, x, y, z)
    wait(device(arch), event)
    return nothing
end

##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    i = @index(Global)
    gi = @index(Group)
    @inbounds ctsdic[x[i], y[i], z[i], gi] += (proc.resp[i] - proc.PS[i]) * ΔT * ac[i]
    @inbounds ctsdoc[x[i], y[i], z[i], gi] += (proc.exu[i] - proc.VDOC[i]) * ΔT * ac[i]
    @inbounds ctsnh4[x[i], y[i], z[i], gi] -= proc.VNH4[i] * ΔT * ac[i]
    @inbounds ctsno3[x[i], y[i], z[i], gi] -= proc.VNO3[i] * ΔT * ac[i]
    @inbounds ctspo4[x[i], y[i], z[i], gi] -= proc.VPO4[i] * ΔT * ac[i]
    # pay attention to ctspop[0,0,0,gi] for non-active individuals
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT, arch::Architecture)
    kernel! = calc_consume_kernel!(device(arch), 1, (size(ac,1),))
    event = kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, proc, ac, x, y, z, ΔT)
    wait(device(arch), event)
    return nothing
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                                   x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global, Linear)
    gi = @index(Group)
    @inbounds ctsdoc[x[i], y[i], z[i], gi] += (plank.Bm[i] + plank.Cq[i]) * lossFracC * plank.ac[i]
    @inbounds ctspoc[x[i], y[i], z[i], gi] += (plank.Bm[i] + plank.Cq[i]) * (1.0 - lossFracC) * plank.ac[i]
    @inbounds ctsdon[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_NC + plank.Nq[i]) * lossFracN * plank.ac[i]
    @inbounds ctspon[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_NC + plank.Nq[i]) * (1.0 - lossFracN) * plank.ac[i]
    @inbounds ctsdop[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_PC + plank.Pq[i]) * lossFracP * plank.ac[i]
    @inbounds ctspop[x[i], y[i], z[i], gi] += (plank.Bm[i] * R_PC + plank.Pq[i]) * (1.0 - lossFracP) * plank.ac[i]
    # pay attention to ctspop[0,0,0,gi] for non-active individuals
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                    x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC, arch::Architecture)
    kernel! = calc_loss_kernel!(device(arch), 1, (size(plank,1),))
    event = kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank,
                    x, y, z, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    wait(device(arch), event)
    return nothing
end

##### calculate PAR field based on Chla field and depth
using KernelAbstractions.Extras.LoopInfo: @unroll
@kernel function calc_par_kernel!(par, chl, PARF, g::Grids, kc, kw)
    @unroll for k in 1:g.Nz
        i, j = @index(Global, NTuple)
        kk = g.Nz - k + 1
        atten = (chl[i,j,kk]/g.V * kc + kw) * g.Δz
        par[i,j,kk] = PARF[i,j] * (1.0 - exp(-atten)) / atten
        PARF[i,j] = PARF[i,j] * exp(-atten)
    end
end
function calc_par!(par, arch::Architecture, chl, PARF, g::Grids, kc, kw)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    event = kernel!(par, chl, PARF, g, kc, kw)
    wait(device(arch), event)
    return nothing
end

