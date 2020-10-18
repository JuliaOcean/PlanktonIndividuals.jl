##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(plank, coord, velos, g::Grids, vel, ΔT)
    find_inds!(plank, coord, plank.ac, g)
    find_vel!(velos, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g, vel.u.data, vel.v.data, vel.w.data)
    find_xᵈ!(velos, plank, g)
    vel1_interpolation!(velos)

    calc_coord!(plank, velos, ΔT)
    periodic_domain!(plank, plank.ac, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(plank, coord, velos, g::Grids, vel₀, vel½, vel₁, ΔT)
    find_inds!(plank, coord, plank.ac, g)
    find_vel!(velos, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g, vel₀.u.data, vel₀.v.data, vel₀.w.data)
    find_xᵈ!(velos, plank, g)
    vel1_interpolation!(velos)

    calc_coord_1!(plank, velos, ΔT)
    periodic_domain!(velos, plank.ac, g)

    find_inds!(velos, coord, plank.ac, g)
    find_vel!(velos, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(velos, velos, g)
    vel2_interpolation!(velos)

    calc_coord_2!(plank, velos, ΔT)
    periodic_domain!(velos, plank.ac, g)

    find_inds!(velos, coord, plank.ac, g)
    find_vel!(velos, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(velos, velos, g)
    vel3_interpolation!(velos)

    calc_coord_3!(plank, velos, ΔT)
    periodic_domain!(velos, plank.ac, g)

    find_inds!(velos, coord, plank.ac, g)
    find_vel!(velos, Int.(coord.x), Int.(coord.y), Int.(coord.z), plank.ac, g, vel₁.u.data, vel₁.v.data, vel₁.w.data)
    find_xᵈ!(velos, velos, g)
    vel4_interpolation!(velos)

    calc_vel_rk4!(velos)
    calc_coord!(plank, velos, ΔT)
    periodic_domain!(plank, plank.ac, g)
end
