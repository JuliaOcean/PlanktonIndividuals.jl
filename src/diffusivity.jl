##### Difference operators with functions
@inline δxᶜᵒᵒ(i, j, k, grid, f::F, args...) where F<:Function = f(i+1, j, k, grid, args...) - f(i, j, k, grid, args...)
@inline δyᵒᶜᵒ(i, j, k, grid, f::F, args...) where F<:Function = f(i, j+1, k, grid, args...) - f(i, j, k, grid, args...)
@inline δzᵒᵒᶜ(i, j, k, grid, f::F, args...) where F<:Function = f(i, j, k+1, grid, args...) - f(i, j, k, grid, args...)

##### Derivative operators
@inline ∂xᶠᵒᵒ(i, j, k, grid, c) = @inbounds (c[i, j, k] - c[i-1, j, k]) / grid.Δx
@inline ∂yᵒᶠᵒ(i, j, k, grid, c) = @inbounds (c[i, j, k] - c[i, j-1, k]) / grid.Δy
@inline ∂zᵒᵒᶠ(i, j, k, grid, c) = @inbounds (c[i, j, k] - c[i, j, k-1]) / grid.Δz

##### Diffusive fluxes
@inline diffusive_flux_x(i, j, k, grid, κˣ::Number, c) = κˣ * grid.Ax * ∂xᶠᵒᵒ(i, j, k, grid, c)
@inline diffusive_flux_y(i, j, k, grid, κʸ::Number, c) = κʸ * grid.Ay * ∂yᵒᶠᵒ(i, j, k, grid, c)
@inline diffusive_flux_z(i, j, k, grid, κᶻ::Number, c) = κᶻ * grid.Az * ∂zᵒᵒᶠ(i, j, k, grid, c)

##### Laplacian diffusion operator
@inline function κ∇²(i, j, k, grid, κˣ, κʸ, κᶻ, c)
    return 1/grid.V * (δxᶜᵒᵒ(i, j, k, grid, diffusive_flux_x, κˣ, c) +
                       δyᵒᶜᵒ(i, j, k, grid, diffusive_flux_y, κʸ, c) +
                       δzᵒᵒᶜ(i, j, k, grid, diffusive_flux_z, κᶻ, c))
end

@inline κ∇²(i, j, k, grid, κ, c) = κ∇²(i, j, k, grid, κ, κ, κ, c)
