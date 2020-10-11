##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(plank, arch::Architecture, g::Grids, vel, ΔT, num::Int64)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel.u.data, vel.v.data, vel.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, 0, num)

    calc_coord!(plank, ΔT, num)
    in_domain!(plank, arch, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(plank, arch::Architecture, g::Grids, vel₀, vel½, vel₁, ΔT, num::Int64)
    find_inds!(plank, arch, g, 12, 0)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₀.u.data, vel₀.v.data, vel₀.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, 0, num)

    calc_intermediate_coord!(plank, ΔT, 0.5, num)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, 1, num)

    calc_intermediate_coord!(plank, ΔT, 0.5, num)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel½.u.data, vel½.v.data, vel½.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, 2, num)

    calc_intermediate_coord!(plank, ΔT, 1.0, num)
    in_domain!(plank, arch, g, 33)

    find_inds!(plank, arch, g, 12, 33)
    find_vel!(plank, Int.(plank[:,13:15]), arch, g, vel₁.u.data, vel₁.v.data, vel₁.w.data)
    find_xᵈ!(plank, Int.(plank[:,13:15]), arch, g)
    vel_interpolation!(plank, 3, num)

    calc_vel_rk4!(plank, num)
    calc_coord!(plank, ΔT, num)
    in_domain!(plank, arch, g)
end
