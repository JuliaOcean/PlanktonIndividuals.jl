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
    plank.x[i] = particle_boundary_condition(plank.x[i], g.xF[g.Hx+1], g.xF[g.Hx+1+g.Nx], TX) * ac[i]
    plank.y[i] = particle_boundary_condition(plank.y[i], g.yF[g.Hy+1], g.yF[g.Hy+1+g.Ny], TY) * ac[i]
    plank.z[i] = particle_boundary_condition(plank.z[i], g.zF[g.Hz+1+g.Nz], g.zF[g.Hz+1], TZ) * ac[i]
end
function particle_boundaries!(plank, ac, g::AbstractGrid, arch::Architecture)
    kernel! = particle_boundaries_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(plank, ac, g)
    wait(device(arch), event)
    return nothing
end

##### calculate intermediate coordinates
@kernel function calc_coord_kernel!(vx, px, u, ac, ΔT, weight::Float64)
    i = @index(Global)
    @inbounds vx[i] = px[i] + weight * u[i] * ΔT * ac[i]
end
function calc_coord!(vx, px, u, ac, ΔT, weight::Float64, arch::Architecture)
    kernel! = calc_coord_kernel!(device(arch), 256, (size(ac,1)))
    event = kernel!(vx, px, u, ac, ΔT, weight)
    wait(device(arch), event)
    return nothing
end

function calc_coord_1!(plank, velos, ΔT, arch::Architecture)
    calc_coord!(velos.x, plank.x, velos.u1, plank.ac, ΔT, 0.5, arch)
    calc_coord!(velos.y, plank.y, velos.v1, plank.ac, ΔT, 0.5, arch)
    calc_coord!(velos.z, plank.z, velos.w1, plank.ac, ΔT, 0.5, arch)
    return nothing
end
function calc_coord_2!(plank, velos, ΔT, arch::Architecture)
    calc_coord!(velos.x, plank.x, velos.u2, plank.ac, ΔT, 0.5, arch)
    calc_coord!(velos.y, plank.y, velos.v2, plank.ac, ΔT, 0.5, arch)
    calc_coord!(velos.z, plank.z, velos.w2, plank.ac, ΔT, 0.5, arch)
    return nothing
end
function calc_coord_3!(plank, velos, ΔT, arch::Architecture)
    calc_coord!(velos.x, plank.x, velos.u3, plank.ac, ΔT, 1.0, arch)
    calc_coord!(velos.y, plank.y, velos.v3, plank.ac, ΔT, 1.0, arch)
    calc_coord!(velos.z, plank.z, velos.w3, plank.ac, ΔT, 1.0, arch)
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

##### calculate final velocities by AB2
@kernel function calc_ab2_kernel!(velos, χ)
    i = @index(Global)
    velos.u1[i] = (1.5 + χ) * velos.u1[i] - (0.5 + χ) * velos.u2[i] 
    velos.v1[i] = (1.5 + χ) * velos.v1[i] - (0.5 + χ) * velos.v2[i] 
    velos.w1[i] = (1.5 + χ) * velos.w1[i] - (0.5 + χ) * velos.w2[i] 
end
function calc_vel_ab2!(velos, χ, arch::Architecture)
    kernel! = calc_ab2_kernel!(device(arch), 256, (size(velos.u1,1)))
    event = kernel!(velos, χ)
    wait(device(arch), event)
    return nothing
end

##### calculate coordinates of each individual
function update_coord!(plank, velos, ΔT, arch::Architecture)
    calc_coord!(plank.x, plank.x, velos.u1, plank.ac, ΔT, 1.0, arch)
    calc_coord!(plank.y, plank.y, velos.v1, plank.ac, ΔT, 1.0, arch)
    calc_coord!(plank.z, plank.z, velos.w1, plank.ac, ΔT, 1.0, arch)
end
