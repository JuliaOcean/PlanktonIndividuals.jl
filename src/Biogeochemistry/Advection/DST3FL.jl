##### halo points ≥ 2 NEEDED!!
const θmax = 1.0e20

##### calculate CFL number: c=uΔt/Lx
@inline CFLx(i, j, k, g::AbstractGrid, u, ΔT) = @inbounds abs(u[i, j, k] * ΔT / ΔxC(i, j, k, g))
@inline CFLy(i, j, k, g::AbstractGrid, v, ΔT) = @inbounds abs(v[i, j, k] * ΔT / ΔyC(i, j, k, g))
@inline CFLz(i, j, k, g::AbstractGrid, w, ΔT) = @inbounds abs(w[i, j, k] * ΔT / ΔzC(i, j, k, g))

##### calculate d₀ and d₁
@inline d0(CFL) = @inbounds (2.0 - CFL) * (1.0 - CFL) / 6.0
@inline d1(CFL) = @inbounds (1.0 - CFL * CFL) / 6.0

##### calculate volume transport, unit: m³/s
@inline Trans_x(i, j, k, g::AbstractGrid, u) = @inbounds Ax(i, j, k, g) * u[i, j, k]
@inline Trans_y(i, j, k, g::AbstractGrid, v) = @inbounds Ay(i, j, k, g) * v[i, j, k]
@inline Trans_z(i, j, k, g::AbstractGrid, w) = @inbounds Az(i, j, k, g) * w[i, j, k]

##### calculate θ⁺ and θ⁻
@inline θx⁺(i, j, k, c, grid) = abs(δx⁰(i, j, k, c, grid))*θmax ≤ abs(δx⁻(i, j, k, c, grid)) ?
    copysign(θmax, δx⁻(i, j, k, c, grid)*δx⁰(i, j, k, c, grid)) : δx⁻(i, j, k, c, grid)/δx⁰(i, j, k, c, grid)

@inline θx⁻(i, j, k, c, grid) = abs(δx⁰(i, j, k, c, grid))*θmax ≤ abs(δx⁺(i, j, k, c, grid)) ?
    copysign(θmax, δx⁺(i, j, k, c, grid)*δx⁰(i, j, k, c, grid)) : δx⁺(i, j, k, c, grid)/δx⁰(i, j, k, c, grid)

@inline θy⁺(i, j, k, c, grid) = abs(δy⁰(i, j, k, c, grid))*θmax ≤ abs(δy⁻(i, j, k, c, grid)) ?
    copysign(θmax, δy⁻(i, j, k, c, grid)*δy⁰(i, j, k, c, grid)) : δy⁻(i, j, k, c, grid)/δy⁰(i, j, k, c, grid)

@inline θy⁻(i, j, k, c, grid) = abs(δy⁰(i, j, k, c, grid))*θmax ≤ abs(δy⁺(i, j, k, c, grid)) ?
    copysign(θmax, δy⁺(i, j, k, c, grid)*δy⁰(i, j, k, c, grid)) : δy⁺(i, j, k, c, grid)/δy⁰(i, j, k, c, grid)

@inline θz⁺(i, j, k, c, grid) = abs(δz⁰(i, j, k, c, grid))*θmax ≤ abs(δz⁻(i, j, k, c, grid)) ?
    copysign(θmax, δz⁻(i, j, k, c, grid)*δz⁰(i, j, k, c, grid)) : δz⁻(i, j, k, c, grid)/δz⁰(i, j, k, c, grid)

@inline θz⁻(i, j, k, c, grid) = abs(δz⁰(i, j, k, c, grid))*θmax ≤ abs(δz⁺(i, j, k, c, grid)) ?
    copysign(θmax, δz⁺(i, j, k, c, grid)*δz⁰(i, j, k, c, grid)) : δz⁺(i, j, k, c, grid)/δz⁰(i, j, k, c, grid)

##### calculate Ψ⁺ and Ψ⁻
@inline Ψx⁺(i, j, k, g::AbstractGrid, u, c, ΔT) = 
            max(0.0, min(min(1.0, d0(CFLx(i, j, k, g, u, ΔT)) + d1(CFLx(i, j, k, g, u, ΔT)) * θx⁺(i, j, k, c, g)), 
                         θx⁺(i, j, k, c, g) * (1.0 - CFLx(i, j, k, g, u, ΔT)) / (CFLx(i, j, k, g, u, ΔT) + 1.0e-20)))

@inline Ψx⁻(i, j, k, g::AbstractGrid, u, c, ΔT) = 
            max(0.0, min(min(1.0, d0(CFLx(i, j, k, g, u, ΔT)) + d1(CFLx(i, j, k, g, u, ΔT)) * θx⁻(i, j, k, c, g)), 
                         θx⁻(i, j, k, c, g) * (1.0 - CFLx(i, j, k, g, u, ΔT)) / (CFLx(i, j, k, g, u, ΔT) + 1.0e-20)))

@inline Ψy⁺(i, j, k, g::AbstractGrid, v, c, ΔT) = 
            max(0.0, min(min(1.0, d0(CFLy(i, j, k, g, v, ΔT)) + d1(CFLy(i, j, k, g, v, ΔT)) * θy⁺(i, j, k, c, g)), 
                         θy⁺(i, j, k, c, g) * (1.0 - CFLy(i, j, k, g, v, ΔT)) / (CFLy(i, j, k, g, v, ΔT) + 1.0e-20)))

@inline Ψy⁻(i, j, k, g::AbstractGrid, v, c, ΔT) =
            max(0.0, min(min(1.0, d0(CFLy(i, j, k, g, v, ΔT)) + d1(CFLy(i, j, k, g, v, ΔT)) * θy⁻(i, j, k, c, g)), 
                         θy⁻(i, j, k, c, g) * (1.0 - CFLy(i, j, k, g, v, ΔT)) / (CFLy(i, j, k, g, v, ΔT) + 1.0e-20)))

@inline Ψz⁺(i, j, k, g::AbstractGrid, w, c, ΔT) = 
            max(0.0, min(min(1.0, d0(CFLz(i, j, k, g, w, ΔT)) + d1(CFLz(i, j, k, g, w, ΔT)) * θz⁺(i, j, k, c, g)), 
                         θz⁺(i, j, k, c, g) * (1.0 - CFLz(i, j, k, g, w, ΔT)) / (CFLz(i, j, k, g, w, ΔT) + 1.0e-20)))

@inline Ψz⁻(i, j, k, g::AbstractGrid, w, c, ΔT) = 
            max(0.0, min(min(1.0, d0(CFLz(i, j, k, g, w, ΔT)) + d1(CFLz(i, j, k, g, w, ΔT)) * θz⁻(i, j, k, c, g)), 
                         θz⁻(i, j, k, c, g) * (1.0 - CFLz(i, j, k, g, w, ΔT)) / (CFLz(i, j, k, g, w, ΔT) + 1.0e-20)))

##### advection flux
@inline adv_flux_x(i, j, k, g::AbstractGrid, u, c, ΔT) = 
    0.5 * (Trans_x(i, j, k, g, u) + abs(Trans_x(i, j, k, g, u))) * (c[i-1, j, k] + Ψx⁺(i, j, k, g, u, c, ΔT) * δx⁰(i, j, k, c, g)) +
    0.5 * (Trans_x(i, j, k, g, u) - abs(Trans_x(i, j, k, g, u))) * (c[i,   j, k] - Ψx⁻(i, j, k, g, u, c, ΔT) * δx⁰(i, j, k, c, g))

@inline adv_flux_y(i, j, k, g::AbstractGrid, v, c, ΔT) =
    0.5 * (Trans_y(i, j, k, g, v) + abs(Trans_y(i, j, k, g, v))) * (c[i, j-1, k] + Ψy⁺(i, j, k, g, v, c, ΔT) * δy⁰(i, j, k, c, g)) +
    0.5 * (Trans_y(i, j, k, g, v) - abs(Trans_y(i, j, k, g, v))) * (c[i, j,   k] - Ψy⁻(i, j, k, g, v, c, ΔT) * δy⁰(i, j, k, c, g))

@inline adv_flux_z(i, j, k, g::AbstractGrid, w, c, ΔT) =
    0.5 * (Trans_z(i, j, k, g, w) + abs(Trans_z(i, j, k, g, w))) * (c[i, j, k  ] + Ψz⁻(i, j, k, g, w, c, ΔT) * δz⁰(i, j, k, c, g)) +
    0.5 * (Trans_z(i, j, k, g, w) - abs(Trans_z(i, j, k, g, w))) * (c[i, j, k-1] - Ψz⁺(i, j, k, g, w, c, ΔT) * δz⁰(i, j, k, c, g))
