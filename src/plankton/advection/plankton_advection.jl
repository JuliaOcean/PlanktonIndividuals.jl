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
function plankton_advectionRK4!(phytos, arch::Architecture, g::Grids, vel₀, vel₁, ΔT)
    vel½ = (u = (vel₀.u .+ vel₁.u) .* 0.5, v = (vel₀.v .+ vel₁.v) .* 0.5, w = (vel₀.w .+ vel₁.w) .* 0.5)
    op_array = adv_op_array_setup(phytos, arch)
    vel_array = vel_array_setup(phytos, arch)
    ind_array = ind_array_setup(phytos, arch)

    find_inds!(ind_array, op_array, arch, g)

    find_vel!(op_array, Int.(ind_array[:,1:3]), arch, g, vel₀.u)
    find_xᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 1)

    find_vel!(op_array, Int.(ind_array[:,4:6]), arch, g, vel₀.v)
    find_yᵈ!(op_array, Int.(ind_array[:,4:6]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 2)

    find_vel!(op_array, Int.(ind_array[:,7:9]), arch, g, vel₀.w)
    find_zᵈ!(op_array, Int.(ind_array[:,7:9]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 3)

    calc_intermediate_coord!(op_array, vel_array, arch, ΔT, 0, 0.5)
    in_domain!(op_array, arch, g, 14)

    update_inds!(ind_array, op_array, arch, g)

    find_vel!(op_array, Int.(ind_array[:,1:3]), arch, g, vel½.u)
    find_xᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 4)

    find_vel!(op_array, Int.(ind_array[:,4:6]), arch, g, vel½.v)
    find_yᵈ!(op_array, Int.(ind_array[:,4:6]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 5)

    find_vel!(op_array, Int.(ind_array[:,7:9]), arch, g, vel½.w)
    find_zᵈ!(op_array, Int.(ind_array[:,7:9]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 6)

    calc_intermediate_coord!(op_array, vel_array, arch, ΔT, 1, 0.5)
    in_domain!(op_array, arch, g, 14)

    update_inds!(ind_array, op_array, arch, g)

    find_vel!(op_array, Int.(ind_array[:,1:3]), arch, g, vel½.u)
    find_xᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 7)

    find_vel!(op_array, Int.(ind_array[:,4:6]), arch, g, vel½.v)
    find_yᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 8)

    find_vel!(op_array, Int.(ind_array[:,7:9]), arch, g, vel½.w)
    find_zᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 9)

    calc_intermediate_coord!(op_array, vel_array, arch, ΔT, 2, 1.0)
    in_domain!(op_array, arch, g, 14)

    update_inds!(ind_array, op_array, arch, g)

    find_vel!(op_array, Int.(ind_array[:,1:3]), arch, g, vel₁.u)
    find_xᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 10)

    find_vel!(op_array, Int.(ind_array[:,4:6]), arch, g, vel₁.v)
    find_yᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 11)

    find_vel!(op_array, Int.(ind_array[:,7:9]), arch, g, vel₁.w)
    find_zᵈ!(op_array, Int.(ind_array[:,1:3]), arch, g)
    vel_interpolation!(vel_array, op_array, arch, g, 12)

    calc_vel_rk4!(vel_array, arch)
    calc_coord!(phytos, vel_array, arch, ΔT)
    in_domain!(phytos, arch, g)
end
