##### Difference operators
@inline δx⁺(i, j, k, c) = @inbounds c[i+1, j, k] - c[i, j, k]
@inline δy⁺(i, j, k, c) = @inbounds c[i, j+1, k] - c[i, j, k]
@inline δz⁺(i, j, k, c) = @inbounds c[i, j, k+1] - c[i, j, k]

@inline δx⁰(i, j, k, c) = @inbounds c[i, j, k] - c[i-1, j, k]
@inline δy⁰(i, j, k, c) = @inbounds c[i, j, k] - c[i, j-1, k]
@inline δz⁰(i, j, k, c) = @inbounds c[i, j, k] - c[i, j, k-1]

##### used only when halo points ≥ 2
@inline δx⁻(i, j, k, c) = @inbounds c[i-1, j, k] - c[i-2, j, k]
@inline δy⁻(i, j, k, c) = @inbounds c[i, j-1, k] - c[i, j-2, k]
@inline δz⁻(i, j, k, c) = @inbounds c[i, j, k-1] - c[i, j, k-2]
#####

##### Difference operators with functions
@inline δx⁺(i, j, k, grid, f::F, args...) where F<:Function = f(i+1, j, k, grid, args...) - f(i, j, k, grid, args...)
@inline δy⁺(i, j, k, grid, f::F, args...) where F<:Function = f(i, j+1, k, grid, args...) - f(i, j, k, grid, args...)
@inline δz⁺(i, j, k, grid, f::F, args...) where F<:Function = f(i, j, k+1, grid, args...) - f(i, j, k, grid, args...)

##### Derivative operators
@inline ∂x⁰(i, j, k, grid, c) = δx⁰(i, j, k, c) / grid.Δx
@inline ∂y⁰(i, j, k, grid, c) = δy⁰(i, j, k, c) / grid.Δy
@inline ∂z⁰(i, j, k, grid, c) = δz⁰(i, j, k, c) / grid.Δz

##### Diffusive fluxes
@inline diffusive_flux_x(i, j, k, grid, κˣ::Number, c) = κˣ * grid.Ax * ∂x⁰(i, j, k, grid, c)
@inline diffusive_flux_y(i, j, k, grid, κʸ::Number, c) = κʸ * grid.Ay * ∂y⁰(i, j, k, grid, c)
@inline diffusive_flux_z(i, j, k, grid, κᶻ::Number, c) = κᶻ * grid.Az * ∂z⁰(i, j, k, grid, c)

##### Laplacian diffusion operator
@inline function κ∇²(i, j, k, grid, κˣ, κʸ, κᶻ, c)
    return 1/grid.V * (δx⁺(i, j, k, grid, diffusive_flux_x, κˣ, c) +
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

function nut_diffusion!(Gcs, arch::Architecture, g, nutrients, κˣ, κʸ, κᶻ, ΔT)
    calc_diffusion_kernel! = calc_diffusion!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    events=[]
    for name in nut_names
        event = calc_diffusion_kernel!(Gcs[name].data, g, κˣ, κʸ, κᶻ, nutrients[name].data, ΔT, dependencies=barrier)
        push!(events, event)
    end

    wait(device(arch), MultiEvent(Tuple(events)))

    return nothing
end
nut_diffusion!(Gcs, arch::Architecture, g, nutrients, κ, ΔT) =
    nut_diffusion!(Gcs, arch::Architecture, g, nutrients, κ, κ, κ, ΔT)
