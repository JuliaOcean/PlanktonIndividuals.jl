##### deal with particles moved out of the domain
@inline function particle_boundary_condition(x, xl, xr, ::Periodic)
    x > xr && return xl + (x  - xr)
    x < xl && return xr - (xl - x) 
    return x
end
@inline function particle_boundary_condition(x, xl, xr, ::Bounded)
    x > xr && return xr - (x - xr) * 1.0f0
    x < xl && return xl + (xl - x) * 1.0f0
    return x
end

@kernel function particle_boundaries_kernel!(particle, ac, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    i = @index(Global)
    @inbounds particle.x[i] = particle_boundary_condition(particle.x[i], 0, g.Nx, TX()) * ac[i]
    @inbounds particle.y[i] = particle_boundary_condition(particle.y[i], 0, g.Ny, TY()) * ac[i]
    @inbounds particle.z[i] = particle_boundary_condition(particle.z[i], 0, g.Nz, TZ()) * ac[i]
end
function particle_boundaries!(particle, ac, g::AbstractGrid, arch::Architecture)
    kernel! = particle_boundaries_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(particle, ac, g)
    return nothing
end

##### find indices (halo points included)
@kernel function find_inds_kernel!(particle, g::AbstractGrid)
    i = @index(Global)
    @inbounds particle.xi[i] = unsafe_trunc(Int, get_xf_index(particle.x[i]) * particle.ac[i]) + g.Hx 
    @inbounds particle.yi[i] = unsafe_trunc(Int, get_yf_index(particle.y[i]) * particle.ac[i]) + g.Hy
    @inbounds particle.zi[i] = unsafe_trunc(Int, get_zf_index(particle.z[i]) * particle.ac[i]) + g.Hz
end
function find_inds!(particle, g::AbstractGrid, arch::Architecture)
    kernel! = find_inds_kernel!(device(arch), 256, (size(particle.ac,1)))
    kernel!(particle, g)
    return nothing
end

##### calculate uvw velocities at (x, y, z)
@kernel function vel_interpolate_kernel!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g::AbstractGrid)
    i = @index(Global)
    @inbounds uₜ[i] = u_itpl(u, x[i], y[i], z[i], ac[i], g) * ac[i]
    @inbounds vₜ[i] = v_itpl(v, x[i], y[i], z[i], ac[i], g) * ac[i]
    @inbounds wₜ[i] = w_itpl(w, x[i], y[i], z[i], ac[i], g) * ac[i]
end

function vel_interpolate!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g::AbstractGrid, arch::Architecture)
    kernel! = vel_interpolate_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_coord_kernel!(velos, particle, u, v, w, ac, ΔT, weight::AbstractFloat)
    i = @index(Global)
    @inbounds velos.x[i] = particle.x[i] + weight * u[i] * ΔT * ac[i]
    @inbounds velos.y[i] = particle.y[i] + weight * v[i] * ΔT * ac[i]
    @inbounds velos.z[i] = particle.z[i] - weight * w[i] * ΔT * ac[i] # index increases when w is negtive
end
function calc_coord!(velos, particle, u, v, w, ac, ΔT, weight::AbstractFloat, arch::Architecture)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(velos, particle, u, v, w, ac, ΔT, weight)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_rk4_kernel!(velos)
    i = @index(Global)
    velos.u1[i] = (velos.u1[i] + 2*velos.u2[i] + 2*velos.u3[i] + velos.u4[i])/6.0f0
    velos.v1[i] = (velos.v1[i] + 2*velos.v2[i] + 2*velos.v3[i] + velos.v4[i])/6.0f0
    velos.w1[i] = (velos.w1[i] + 2*velos.w2[i] + 2*velos.w3[i] + velos.w4[i])/6.0f0
end
function calc_vel_rk4!(velos, arch::Architecture)
    kernel! = calc_rk4_kernel!(device(arch), 256, (size(velos.u1,1)))
    kernel!(velos)
    return nothing
end

##### update coordinates of each individual using RK4 integration
function particle_advection!(particle, velos, g::AbstractGrid, vel₀, vel½, vel₁, ΔT, arch::Architecture)
    vel_interpolate!(velos.u1, velos.v1, velos.w1, particle.x, particle.y, particle.z, particle.ac, 
                     vel₀.u.data, vel₀.v.data, vel₀.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .= velos.u1
    velos.v2 .= velos.v1
    velos.w2 .= velos.w1

    calc_coord!(velos, particle, velos.u1, velos.v1, velos.w1, particle.ac, ΔT, 0.5f0, arch)
    particle_boundaries!(velos, particle.ac, g, arch)

    ##### stage 2
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, particle.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1 .* 2
    velos.v2 .+= velos.v1 .* 2
    velos.w2 .+= velos.w1 .* 2

    calc_coord!(velos, particle, velos.u1, velos.v1, velos.w1, particle.ac, ΔT, 0.5f0, arch)
    particle_boundaries!(velos, particle.ac, g, arch)

    ##### stage 3
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, particle.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1 .* 2
    velos.v2 .+= velos.v1 .* 2
    velos.w2 .+= velos.w1 .* 2

    calc_coord!(velos, particle, velos.u1, velos.v1, velos.w1, particle.ac, ΔT, 1.0f0, arch)
    particle_boundaries!(velos, particle.ac, g, arch)

    ##### stage 4
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, particle.ac, 
                     vel₁.u.data, vel₁.v.data, vel₁.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1
    velos.v2 .+= velos.v1
    velos.w2 .+= velos.w1

    ##### calculate final velocities
    velos.u2 .= velos.u2 ./ 6
    velos.v2 .= velos.v2 ./ 6
    velos.w2 .= velos.w2 ./ 6

    calc_coord!(particle, particle, velos.u2, velos.v2, velos.w2, particle.ac, ΔT, 1.0f0, arch)
    particle_boundaries!(particle, particle.ac, g, arch)
    
    return nothing
end