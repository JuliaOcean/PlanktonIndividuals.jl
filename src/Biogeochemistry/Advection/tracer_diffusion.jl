##### Diffusive fluxes
@inline diffusive_flux_x(i, j, k, grid, κˣ::Number, c) = κˣ * Ax(i, j, k, grid) * ∂x⁰(i, j, k, grid, c)
@inline diffusive_flux_y(i, j, k, grid, κʸ::Number, c) = κʸ * Ay(i, j, k, grid) * ∂y⁰(i, j, k, grid, c)
@inline diffusive_flux_z(i, j, k, grid, κᶻ::Number, c) = κᶻ * Az(i, j, k, grid) * ∂z⁰(i, j, k, grid, c)

##### Laplacian diffusion operator
@inline function κ∇²(i, j, k, grid, κˣ, κʸ, κᶻ, c)
    return 1.0f0/volume(i, j, k, grid) * (δx⁺(i, j, k, grid, diffusive_flux_x, κˣ, c) +
                                          δy⁺(i, j, k, grid, diffusive_flux_y, κʸ, c) +
                                          δz⁺(i, j, k, grid, diffusive_flux_z, κᶻ, c))
end

##### calculate the tendency by diffusion for tracer c
@kernel function calc_diffusion!(Gc, grid, κˣ, κʸ, κᶻ, c, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + grid.Hx
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[ii, jj, kk] = Gc[ii, jj, kk] + κ∇²(ii, jj, kk, grid, κˣ, κʸ, κᶻ, c) * ΔT
end

function tracer_diffusion!(Gcs, arch::Architecture, g, tracers, κˣ, κʸ, κᶻ, ΔT)
    calc_diffusion_kernel! = calc_diffusion!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        calc_diffusion_kernel!(Gcs[name].data, g, κˣ, κʸ, κᶻ, tracers[name].data, ΔT)
    end

    return nothing
end
