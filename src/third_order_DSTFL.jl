###################################################################################
# Compute advective flux of tracers using 3rd order DST Scheme with flux limiting #
# Due to free surface, tracers will not be conserved at each time step.           #
# Surface flux(τⁿ[:,:,1]*Az*wFld[:,:,1]) need to be recorded to check tracer      #
# conservation.                                                                   #
# There will still be some tiny negative values in tracer field because of multi- #
# dimensional advection.                                                          #
# velocities should not be too high, otherwise big negative values will occur     #
###################################################################################

##### halo points ≥ 2 NEEDED!!

const θmax = 1.0e20

##### calculate CFL number: c=uΔt/Lx
@inline CFLx(i, j, k, g::grids, u, ΔT) = @inbounds abs(u[i, j, k] * ΔT / g.Δx)
@inline CFLy(i, j, k, g::grids, v, ΔT) = @inbounds abs(v[i, j, k] * ΔT / g.Δy)
@inline CFLz(i, j, k, g::grids, w, ΔT) = @inbounds abs(w[i, j, k] * ΔT / g.Δz)

##### calculate d₀ and d₁
@inline d0(CFL) = @inbounds (2.0 - CFL) * (1.0 - CFL) / 6.0
@inline d1(CFL) = @inbounds (1.0 - CFL * CFL) / 6.0

##### calculate volume transport, unit: m³/s
@inline Trans_x(i, j, k, g::grids, u) = @inbounds g.Ax * u[i, j, k]
@inline Trans_y(i, j, k, g::grids, v) = @inbounds g.Ay * v[i, j, k]
@inline Trans_z(i, j, k, g::grids, w) = @inbounds g.Az * w[i, j, k]

##### calculate θ⁺ and θ⁻
@inline θx⁺(i, j, k, q) = abs(δx⁰(i, j, k, q))*θmax ≤ abs(δx⁻(i, j, k, q)) ?
    copysign(θmax, δx⁻(i, j, k, q)*δx⁰(i, j, k, q)) : δx⁻(i, j, k, q)/δx⁰(i, j, k, q)

@inline θx⁻(i, j, k, q) = abs(δx⁰(i, j, k, q))*θmax ≤ abs(δx⁺(i, j, k, q)) ?
    copysign(θmax, δx⁺(i, j, k, q)*δx⁰(i, j, k, q)) : δx⁺(i, j, k, q)/δx⁰(i, j, k, q)

@inline θy⁺(i, j, k, q) = abs(δy⁰(i, j, k, q))*θmax ≤ abs(δy⁻(i, j, k, q)) ?
    copysign(θmax, δy⁻(i, j, k, q)*δy⁰(i, j, k, q)) : δy⁻(i, j, k, q)/δy⁰(i, j, k, q)

@inline θy⁻(i, j, k, q) = abs(δy⁰(i, j, k, q))*θmax ≤ abs(δy⁺(i, j, k, q)) ?
    copysign(θmax, δy⁺(i, j, k, q)*δy⁰(i, j, k, q)) : δy⁺(i, j, k, q)/δy⁰(i, j, k, q)

@inline θz⁺(i, j, k, q) = abs(δz⁰(i, j, k, q))*θmax ≤ abs(δz⁻(i, j, k, q)) ?
    copysign(θmax, δz⁻(i, j, k, q)*δz⁰(i, j, k, q)) : δz⁻(i, j, k, q)/δz⁰(i, j, k, q)

@inline θz⁻(i, j, k, q) = abs(δz⁰(i, j, k, q))*θmax ≤ abs(δz⁺(i, j, k, q)) ?
    copysign(θmax, δz⁺(i, j, k, q)*δz⁰(i, j, k, q)) : δz⁺(i, j, k, q)/δz⁰(i, j, k, q)

##### calculate Ψ⁺ and Ψ⁻
@inline Ψx⁺(i, j, k, g::grids, u, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLx(i, j, k, g::grids, u, ΔT)) +
                                                            d1(CFLx(i, j, k, g::grids, u, ΔT)) *
                                                            θx⁺(i, j, k, q)
                                                            ),
                                                        θx⁺(i, j, k, q) * (1.0 - CFLx(i, j, k, g::grids, u, ΔT)) /
                                                        (CFLx(i, j, k, g::grids, u, ΔT) + 1.0e-20)
                                                        )
                                               )
@inline Ψx⁻(i, j, k, g::grids, u, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLx(i, j, k, g::grids, u, ΔT)) +
                                                            d1(CFLx(i, j, k, g::grids, u, ΔT)) *
                                                            θx⁻(i, j, k, q)
                                                            ),
                                                        θx⁻(i, j, k, q) * (1.0 - CFLx(i, j, k, g::grids, u, ΔT)) /
                                                        (CFLx(i, j, k, g::grids, u, ΔT) + 1.0e-20)
                                                        )
                                               )

@inline Ψy⁺(i, j, k, g::grids, v, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLy(i, j, k, g::grids, v, ΔT)) +
                                                            d1(CFLy(i, j, k, g::grids, v, ΔT)) *
                                                            θy⁺(i, j, k, q)
                                                            ),
                                                        θy⁺(i, j, k, q) * (1.0 - CFLy(i, j, k, g::grids, v, ΔT)) /
                                                        (CFLy(i, j, k, g::grids, v, ΔT) + 1.0e-20)
                                                        )
                                               )
@inline Ψy⁻(i, j, k, g::grids, v, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLy(i, j, k, g::grids, v, ΔT)) +
                                                            d1(CFLy(i, j, k, g::grids, v, ΔT)) *
                                                            θy⁻(i, j, k, q)
                                                            ),
                                                        θy⁻(i, j, k, q) * (1.0 - CFLy(i, j, k, g::grids, v, ΔT)) /
                                                        (CFLy(i, j, k, g::grids, v, ΔT) + 1.0e-20)
                                                        )
                                               )

@inline Ψz⁺(i, j, k, g::grids, w, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLz(i, j, k, g::grids, w, ΔT)) +
                                                            d1(CFLz(i, j, k, g::grids, w, ΔT)) *
                                                            θz⁺(i, j, k, q)
                                                            ),
                                                        θz⁺(i, j, k, q) * (1.0 - CFLz(i, j, k, g::grids, w, ΔT)) /
                                                        (CFLz(i, j, k, g::grids, w, ΔT) + 1.0e-20)
                                                        )
                                               )
@inline Ψz⁻(i, j, k, g::grids, w, q, ΔT) = max(0.0, min(min(1.0,
                                                            d0(CFLz(i, j, k, g::grids, w, ΔT)) +
                                                            d1(CFLz(i, j, k, g::grids, w, ΔT)) *
                                                            θz⁻(i, j, k, q)
                                                            ),
                                                        θz⁻(i, j, k, q) * (1.0 - CFLz(i, j, k, g::grids, w, ΔT)) /
                                                        (CFLz(i, j, k, g::grids, w, ΔT) + 1.0e-20)
                                                        )
                                               )

##### advection flux
@inline adv_flux_x(i, j, k, g::grids, u, q, ΔT) =
    0.5 * (Trans_x(i, j, k, g::grids, u) + abs(Trans_x(i, j, k, g::grids, u))) *
    (q[i-1, j, k] + Ψx⁺(i, j, k, g::grids, u, q, ΔT) * δx⁰(i, j, k, q)) +
    0.5 * (Trans_x(i, j, k, g::grids, u) - abs(Trans_x(i, j, k, g::grids, u))) *
    (q[i,   j, k] + Ψx⁻(i, j, k, g::grids, u, q, ΔT) * δx⁰(i, j, k, q))

@inline adv_flux_y(i, j, k, g::grids, v, q, ΔT) =
    0.5 * (Trans_y(i, j, k, g::grids, v) + abs(Trans_y(i, j, k, g::grids, v))) *
    (q[i, j-1, k] + Ψy⁺(i, j, k, g::grids, v, q, ΔT) * δy⁰(i, j, k, q)) +
    0.5 * (Trans_y(i, j, k, g::grids, v) - abs(Trans_y(i, j, k, g::grids, v))) *
    (q[i, j,   k] + Ψy⁻(i, j, k, g::grids, v, q, ΔT) * δy⁰(i, j, k, q))

@inline adv_flux_z(i, j, k, g::grids, w, q, ΔT) =
    0.5 * (Trans_z(i, j, k, g::grids, w) + abs(Trans_z(i, j, k, g::grids, w))) *
    (q[i, j, k-1] + Ψz⁺(i, j, k, g::grids, w, q, ΔT) * δz⁰(i, j, k, q)) +
    0.5 * (Trans_z(i, j, k, g::grids, w) - abs(Trans_z(i, j, k, g::grids, w))) *
    (q[i, j,   k] + Ψz⁻(i, j, k, g::grids, w, q, ΔT) * δw⁰(i, j, k, q))
