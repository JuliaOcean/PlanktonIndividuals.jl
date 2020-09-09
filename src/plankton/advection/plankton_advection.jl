##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(phytos, arch::Architecture, g::Grids, vel, ΔT)
    op_array = adv_op_array_setup(phytos, arch)
    vel_array = vel_array_setup(phytos, arch)

    find_inds!(op_array, arch, g)

    find_u_vels!(op_array, arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 1)

    find_v_vels!(op_array, arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 2)

    find_w_vels!(op_array, arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 3)

    calc_coord!(phytos, vel_array, arch, g, ΔT)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(phytos, arch::Architecture, g::Grids, vel₀, vel₁, ΔT)
    vel½ = (u = (vel₀.u .+ vel₁.u) .* 0.5, v = (vel₀.v .+ vel₁.v) .* 0.5, w = (vel₀.w .+ vel₁.w) .* 0.5)
    op_array = adv_op_array_setup(phytos, arch)
    vel_array = vel_array_setup(phytos, arch)

    find_inds!(op_array, arch, g)
    find_u_vels!(op_array, arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 1)
    find_v_vels!(op_array, arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 2)
    find_w_vels!(op_array, arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 3)
    calc_intermediate_coord!(op_array, vel_array, arch, g, ΔT, 0)

    update_inds!(op_array, arch, g)
    find_u_vels!(op_array, arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 4)
    find_v_vels!(op_array, arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 5)
    find_w_vels!(op_array, arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 6)
    calc_intermediate_coord!(op_array, vel_array, arch, g, ΔT, 1)

    update_inds!(op_array, arch, g)
    find_u_vels!(op_array, arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 7)
    find_v_vels!(op_array, arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 8)
    find_w_vels!(op_array, arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 9)
    calc_intermediate_coord!(op_array, vel_array, arch, g, ΔT, 2)

    update_inds!(op_array, arch, g)
    find_u_vels!(op_array, arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 10)
    find_v_vels!(op_array, arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 11)
    find_w_vels!(op_array, arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 12)

    calc_vel_rk4!(vel_array, arch)
    calc_coord!(phytos, vel_array, arch, g, ΔT)
end
