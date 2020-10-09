using KernelAbstractions.Extras.LoopInfo: @unroll
##### calculate Chla field based on the status of plankton individuals
@kernel function acc_chla_field_kernel!(chl, pop, plank, inds::AbstractArray{Int64,2})
    @unroll for i in 1:size(plank,1)
        if plank[i,58] == 1.0
            @inbounds xi = inds[i,1]
            @inbounds yi = inds[i,2]
            @inbounds zi = inds[i,3]
            @inbounds chl[xi, yi, zi] += plank[i,10]
            @inbounds pop[xi, yi, zi] += 1.0
        end
    end
end
function acc_chla_field!(chl, pop, plank, inds::AbstractArray{Int64,2}, arch::Architecture)
    kernel! = acc_chla_field_kernel!(device(arch), 256, (1,))
    event = kernel!(chl, pop, plank, inds)
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

