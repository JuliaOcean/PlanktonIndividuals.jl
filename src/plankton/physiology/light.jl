using KernelAbstractions.Extras.LoopInfo: @unroll
##### calculate Chla field based on the status of plankton individuals
@kernel function calc_chla_field_kernel!(chl, phytos, g::Grids)
    i = @index(Global, Linear)
    xi = find_xF_ind(phytos[i,1], g) |> Int
    yi = find_yF_ind(phytos[i,2], g) |> Int
    zi = find_zF_ind(phytos[i,3], g) |> Int
    chl[xi, yi, zi] = chl[xi, yi, zi] + phytos[i,9]
end
function calc_chla_field!(chl, arch::Architecture, phytos, g::Grids)
    kernel! = calc_chla_field_kernel!(device(arch), 256, (size(phytos,1),))
    event = kernel!(chl, phytos, g)
    wait(device(arch), event)
    chl .= chl ./ g.V
    return nothing
end

##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, chl, g::Grids, surf_par, kc, kw)
    i, j = @index(Global, NTuple)
    @unroll for k in grid.Nz:-1:1
        par[i,j,k] = surf_par[i,j] * (1.0 - exp(-(chl[i,j,k] * kc + kw) * g.Δz)) / ((chl[i,j,k] * kc + kw) * g.Δz)
        surf_par[i,j] = surf_par[i,j] * exp(-(chl[i,j,k] * kc + kw) * g.Δz)
    end
end
function calc_par!(par, arch::Architecture, chl, g::Grids, surf_par, kc, kw)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    event = kernel!(par, chl, g, surf_par, kc, kw)
    wait(device(arch), event)
    return nothing
end

