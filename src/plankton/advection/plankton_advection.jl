##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(plank, arch::Architecture, g::Grids, vel, ΔT)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel.u.data, vel.v.data, vel.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 0)

    calc_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(plank, arch::Architecture, g::Grids, vel₀, vel½, vel₁, ΔT)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₀.u.data, vel₀.v.data, vel₀.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, 0)

    calc_1st_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, 1)

    calc_2nd_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, 2)

    calc_3rd_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₁.u.data, vel₁.v.data, vel₁.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, 3)

    calc_vel_rk4!(plank, arch)
    calc_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g)
end
