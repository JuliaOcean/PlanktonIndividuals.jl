##### Difference operators
@inline δx⁺(i, j, k, c, grid) = @inbounds (c[i+1, j, k] - c[i, j, k]) * grid.landmask[i+1, j,   k] * grid.landmask[i, j, k]
@inline δy⁺(i, j, k, c, grid) = @inbounds (c[i, j+1, k] - c[i, j, k]) * grid.landmask[i,   j+1, k] * grid.landmask[i, j, k]
@inline δz⁺(i, j, k, c, grid) = @inbounds (c[i, j, k] - c[i, j, k+1]) * grid.landmask[i,   j, k+1] * grid.landmask[i, j, k]

@inline δx⁰(i, j, k, c, grid) = @inbounds (c[i, j, k] - c[i-1, j, k]) * grid.landmask[i, j, k] * grid.landmask[i-1, j,   k]
@inline δy⁰(i, j, k, c, grid) = @inbounds (c[i, j, k] - c[i, j-1, k]) * grid.landmask[i, j, k] * grid.landmask[i,   j-1, k]
@inline δz⁰(i, j, k, c, grid) = @inbounds (c[i, j, k-1] - c[i, j, k]) * grid.landmask[i, j, k] * grid.landmask[i,   j, k-1]

##### used only when halo points ≥ 2
@inline δx⁻(i, j, k, c, grid) = @inbounds (c[i-1, j, k] - c[i-2, j, k]) * grid.landmask[i-1, j,   k] * grid.landmask[i-2, j,   k]
@inline δy⁻(i, j, k, c, grid) = @inbounds (c[i, j-1, k] - c[i, j-2, k]) * grid.landmask[i,   j-1, k] * grid.landmask[i,   j-2, k]
@inline δz⁻(i, j, k, c, grid) = @inbounds (c[i, j, k-2] - c[i, j, k-1]) * grid.landmask[i,   j, k-1] * grid.landmask[i,   j, k-2]
#####

##### Difference operators with functions
@inline δx⁺(i, j, k, grid, f::F, args...) where F<:Function = (f(i+1, j, k, grid, args...) - f(i, j, k, grid, args...)) *
                                                                grid.landmask[i+1, j,   k] * grid.landmask[i, j, k]
@inline δy⁺(i, j, k, grid, f::F, args...) where F<:Function = (f(i, j+1, k, grid, args...) - f(i, j, k, grid, args...)) *
                                                                grid.landmask[i,   j+1, k] * grid.landmask[i, j, k]
@inline δz⁺(i, j, k, grid, f::F, args...) where F<:Function = (f(i, j, k, grid, args...) - f(i, j, k+1, grid, args...)) *
                                                                grid.landmask[i,   j, k+1] * grid.landmask[i, j, k]

##### Derivative operators
@inline ∂x⁰(i, j, k, grid, c) = δx⁰(i, j, k, c, grid) / ΔxF(i, j, k, grid)
@inline ∂y⁰(i, j, k, grid, c) = δy⁰(i, j, k, c, grid) / ΔyF(i, j, k, grid)
@inline ∂z⁰(i, j, k, grid, c) = δz⁰(i, j, k, c, grid) / ΔzF(i, j, k, grid)