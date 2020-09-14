##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(plank, arch::Architecture, g::Grids, vel, ΔT)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel.u, vel.v, vel.w)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 0)

    calc_coord!(plank, plank, arch, ΔT)
    in_domain!(plank, arch, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(plank, arch::Architecture, g::Grids, vel₀, vel½, vel₁, ΔT)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₀.u, vel₀.v, vel₀.w)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 0)

    calc_1st_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 1)

    calc_2nd_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 2)

    calc_3rd_intermediate_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₁.u, vel₁.v, vel₁.w)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, arch, g, 3)

    calc_vel_rk4!(plank, arch)
    calc_coord!(plank, arch, ΔT)
    in_domain!(plank, arch, g)
end
