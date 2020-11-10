using KernelAbstractions.Extras.LoopInfo: @unroll
##### calculate Chla and individual counts based on the status of plankton individuals
@kernel function acc_counts_kernel!(ctschl, ctsppl, chl, ac, x, y, z)
    i = @index(Global)
    gi = @index(Group)
    @inbounds ctschl[x[i], y[i], z[i], gi] += copy(chl[i]) * ac[i]
    @inbounds ctsppl[x[i], y[i], z[i], gi] += 1.0 * ac[i]
    # pay attention to cumpop[0,0,0,gi] for non-active individuals
end
function acc_counts!(ctschl, ctsppl, chl, ac, x, y, z, arch::Architecture)
    kernel! = acc_counts_kernel!(device(arch), 1, (size(ac,1),))
    event = kernel!(ctschl, ctsppl, chl, ac, x, y, z)
    wait(device(arch), event)
    return nothing
end

##### calculate PAR field based on Chla field and depth
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

