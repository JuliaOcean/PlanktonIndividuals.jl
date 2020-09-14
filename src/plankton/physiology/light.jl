using KernelAbstractions.Extras.LoopInfo: @unroll
##### calculate Chla field based on the status of plankton individuals
@kernel function acc_chla_field_kernel!(chl, plank, inds::AbstractArray{Int64,2})
    i = @index(Global, Linear)
    xi = inds[i,1]
    yi = inds[i,2]
    zi = inds[i,3]
    chl[xi, yi, zi] = chl[xi, yi, zi] + plank[i,10]
end
function acc_chla_field!(chl, plank, inds::AbstractArray{Int64,2}, arch::Architecture)
    kernel! = acc_chla_field_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(chl, plank, inds)
    wait(device(arch), event)
    return nothing
end

##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, chl, g::Grids, surf_par, kc, kw)
    @unroll for k in 1:g.Nz
        i, j = @index(Global, NTuple)
        kk = g.Nz - k + 1
        atten = exp(-(chl[i,j,kk]/g.V * kc + kw) * g.Δz)
        par[i,j,kk] = surf_par[i,j] * (1.0 - atten) / ((chl[i,j,kk]/g.V * kc + kw) * g.Δz)
        surf_par[i,j] = surf_par[i,j] * atten
    end

end
function calc_par!(par, arch::Architecture, chl, g::Grids, surf_par, kc, kw)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    event = kernel!(par, chl, g, surf_par, kc, kw)
    wait(device(arch), event)
    return nothing
end

