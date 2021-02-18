##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(plank, velos, g::Grids, vel, ΔT, arch::Architecture)
    vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
                     vel.u.data, vel.v.data, vel.w.data, g, arch)

    update_coord!(plank, velos, ΔT, arch)
    periodic_domain!(plank, plank.ac, g, arch)
end

##### update coordinates of each individual using RK4 integration
function plankton_advection!(plank, velos, g::Grids, vel₀, vel½, vel₁, ΔT, arch::Architecture)
    vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
                     vel₀.u.data, vel₀.v.data, vel₀.w.data, g, arch)
    calc_coord_1!(plank, velos, ΔT, arch)
    periodic_domain!(velos, plank.ac, g, arch)

    ##### stage 2
    vel_interpolate!(velos.u2, velos.v2, velos.w2, velos.x, velos.y, velos.z, plank.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)
    calc_coord_2!(plank, velos, ΔT, arch)
    periodic_domain!(velos, plank.ac, g, arch)

    ##### stage 3
    vel_interpolate!(velos.u3, velos.v3, velos.w3, velos.x, velos.y, velos.z, plank.ac, 
                     vel½.u.data, vel½.v.data, vel½.w.data, g, arch)
    calc_coord_3!(plank, velos, ΔT, arch)
    periodic_domain!(velos, plank.ac, g, arch)

    ##### stage 4
    vel_interpolate!(velos.u4, velos.v4, velos.w4, velos.x, velos.y, velos.z, plank.ac, 
                     vel₁.u.data, vel₁.v.data, vel₁.w.data, g, arch)

    calc_vel_rk4!(velos, arch)
    update_coord!(plank, velos, ΔT, arch)
    periodic_domain!(plank, plank.ac, g, arch)
end

##### update coordinates of each individual using quasiAB2 integration
function plankton_advection!(plank, velos, g::Grids, χ, vel₁, ΔT, arch::Architecture)
    ##### time step n
    vel_interpolate!(velos.u1, velos.v1, velos.w1, plank.x, plank.y, plank.z, plank.ac, 
                     vel₁.u.data, vel₁.v.data, vel₁.w.data, g, arch)

    ##### AB2
    calc_vel_ab2!(velos, χ, arch)
    update_coord!(plank, velos, ΔT, arch)

    ##### time step n-1
    velos.u2 .= velos.u1
    velos.v2 .= velos.v1
    velos.w2 .= velos.w1

    return nothing
end