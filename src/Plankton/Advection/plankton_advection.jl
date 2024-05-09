##### update coordinates of each individual using RK4 integration
function plankton_advection!(plank, velos, g::AbstractGrid, vel₀, vel½, vel₁, ΔT, arch::Architecture)
    vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
                     vel₀.u.data, vel₀.v.data, vel₀.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .= velos.u1
    velos.v2 .= velos.v1
    velos.w2 .= velos.w1

    calc_coord!(velos, plank, velos.u1, velos.v1, velos.w1, plank.ac, ΔT, 0.5f0, arch)
    particle_boundaries!(velos, plank.ac, g, arch)

    ##### stage 2
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, plank.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1 .* 2
    velos.v2 .+= velos.v1 .* 2
    velos.w2 .+= velos.w1 .* 2

    calc_coord!(velos, plank, velos.u1, velos.v1, velos.w1, plank.ac, ΔT, 0.5f0, arch)
    particle_boundaries!(velos, plank.ac, g, arch)

    ##### stage 3
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, plank.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1 .* 2
    velos.v2 .+= velos.v1 .* 2
    velos.w2 .+= velos.w1 .* 2

    calc_coord!(velos, plank, velos.u1, velos.v1, velos.w1, plank.ac, ΔT, 1.0f0, arch)
    particle_boundaries!(velos, plank.ac, g, arch)

    ##### stage 4
    vel_interpolate!(velos.u1, velos.v1, velos.w1, velos.x, velos.y, velos.z, plank.ac, 
                     vel₁.u.data, vel₁.v.data, vel₁.w.data, g, arch)

    ##### add up intermediate velocities
    velos.u2 .+= velos.u1
    velos.v2 .+= velos.v1
    velos.w2 .+= velos.w1

    ##### calculate final velocities
    velos.u2 .= velos.u2 ./ 6
    velos.v2 .= velos.v2 ./ 6
    velos.w2 .= velos.w2 ./ 6

    calc_coord!(plank, plank, velos.u2, velos.v2, velos.w2, plank.ac, ΔT, 1.0f0, arch)
    particle_boundaries!(plank, plank.ac, g, arch)
    
    return nothing
end

# ##### update coordinates of each individual using quasiAB2 integration
# function plankton_advection!(plank, velos, g::AbstractGrid, χ, vel₁, ΔT, arch::Architecture)
#     ##### time step n
#     vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
#                      vel₁.u.data, vel₁.v.data, vel₁.w.data, g, arch)

#     ##### AB2
#     calc_vel_ab2!(velos, χ, arch)
#     calc_coord!(plank, plank, velos.u1, velos.v1, velos.w1, plank.ac, ΔT, 1.0, arch)

#     ##### time step n-1
#     velos.u2 .= velos.u1
#     velos.v2 .= velos.v1
#     velos.w2 .= velos.w1

#     return nothing
# end

# ##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
# function plankton_advection!(plank, velos, g::AbstractGrid, vel, ΔT, arch::Architecture)
#     vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
#                      vel.u.data, vel.v.data, vel.w.data, g, arch)

#     calc_coord!(plank, plank, velos.u1, velos.v1, velos.w1, plank.ac, ΔT, 1.0, arch)
#     particle_boundaries!(plank, plank.ac, g, arch)

#     return nothing
# end
