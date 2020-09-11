##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(phytos, arch::Architecture, g::Grids, vel, ΔT)
    op_array = adv_op_array_setup(phytos, arch)
    vel_array = vel_array_setup(phytos, arch)
    ind_array = ind_array_setup(phytos, arch)

    find_inds!(ind_array, op_array, arch, g)

    find_vel!(op_array, ind_array[:,1:3], arch, g, vel.u)
    vel_interpolation!(vel_array, op_array, arch, g, 1)

    find_vel!(op_array, ind_array[:,4:6], arch, g, vel.v)
    vel_interpolation!(vel_array, op_array, arch, g, 2)

    find_vel!(op_array, ind_array[:,7:9], arch, g, vel.w)
    vel_interpolation!(vel_array, op_array, arch, g, 3)

    calc_coord!(phytos, vel_array, arch, g, ΔT)
    in_domain!(phytos, arch, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(phytos, ope, arch::Architecture, g::Grids, vel₀, vel½, vel₁, ΔT)
    update_op_array!(phytos, ope, arch)

    find_inds!(ope, arch, g)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel₀.u, vel₀.v, vel₀.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 0)

    calc_1st_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    update_inds!(ope, arch, g)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 1)

    calc_2nd_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    update_inds!(ope, arch, g)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 2)

    calc_3rd_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    update_inds!(ope, arch, g)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel₁.u, vel₁.v, vel₁.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 3)

    calc_vel_rk4!(ope, arch)
    calc_coord!(phytos, ope, arch, ΔT)
    in_domain!(phytos, arch, g)
end
