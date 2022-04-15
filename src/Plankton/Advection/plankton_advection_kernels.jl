##### deal with particles moved out of the domain
@inline function particle_boundary_condition(x, xl, xr, ::Periodic)
    x > xr && return xl + (x  - xr)
    x < xl && return xr - (xl - x) 
    return x
end
@inline function particle_boundary_condition(x, xl, xr, ::Bounded)
    x > xr && return xr - (x - xr) * 1.0
    x < xl && return xl + (xl - x) * 1.0
    return x
end

@kernel function particle_boundaries_kernel!(plank, ac, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    i = @index(Global)
    @inbounds plank.x[i] = particle_boundary_condition(plank.x[i], 0, g.Nx, TX()) * ac[i]
    @inbounds plank.y[i] = particle_boundary_condition(plank.y[i], 0, g.Ny, TY()) * ac[i]
    @inbounds plank.z[i] = particle_boundary_condition(plank.z[i], 0, g.Nz, TZ()) * ac[i]
end
function particle_boundaries!(plank, ac, g::AbstractGrid, arch::Architecture)
    kernel! = particle_boundaries_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(plank, ac, g)
    wait(device(arch), event)
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
    event = kernel!(uₜ, vₜ, wₜ, x, y, z, ac, u, v, w, g)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_coord_kernel!(velos, plank, u, v, w, ac, ΔT, weight::Float64)
    i = @index(Global)
    @inbounds velos.x[i] = plank.x[i] + weight * u[i] * ΔT * ac[i]
    @inbounds velos.y[i] = plank.y[i] + weight * v[i] * ΔT * ac[i]
    @inbounds velos.z[i] = plank.z[i] - weight * w[i] * ΔT * ac[i] # index increases when w is negtive
end
function calc_coord!(velos, plank, u, v, w, ac, ΔT, weight::Float64, arch::Architecture)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(velos, plank, u, v, w, ac, ΔT, weight)
    wait(device(arch), event)
    return nothing
end

##### calculate final velocities by RK4
@kernel function calc_rk4_kernel!(velos)
    i = @index(Global)
    velos.u1[i] = (velos.u1[i] + 2*velos.u2[i] + 2*velos.u3[i] + velos.u4[i])/6 
    velos.v1[i] = (velos.v1[i] + 2*velos.v2[i] + 2*velos.v3[i] + velos.v4[i])/6 
    velos.w1[i] = (velos.w1[i] + 2*velos.w2[i] + 2*velos.w3[i] + velos.w4[i])/6 
end
function calc_vel_rk4!(velos, arch::Architecture)
    kernel! = calc_rk4_kernel!(device(arch), 256, (size(velos.u1,1)))
    event = kernel!(velos)
    wait(device(arch), event)
    return nothing
end

# ##### calculate final velocities by AB2
# @kernel function calc_ab2_kernel!(velos, χ)
#     i = @index(Global)
#     velos.u1[i] = (1.5 + χ) * velos.u1[i] - (0.5 + χ) * velos.u2[i] 
#     velos.v1[i] = (1.5 + χ) * velos.v1[i] - (0.5 + χ) * velos.v2[i] 
#     velos.w1[i] = (1.5 + χ) * velos.w1[i] - (0.5 + χ) * velos.w2[i] 
# end
# function calc_vel_ab2!(velos, χ, arch::Architecture)
#     kernel! = calc_ab2_kernel!(device(arch), 256, (size(velos.u1,1)))
#     event = kernel!(velos, χ)
#     wait(device(arch), event)
#     return nothing
# end