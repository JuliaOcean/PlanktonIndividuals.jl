##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(phytos, arch::Architecture, g::Grids, vel, ΔT)
    update_adv_ope!(phytos, ope, arch)

    find_inds!(ope, arch, g, 1, 0)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel.u, vel.v, vel.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 0)

    calc_coord!(phytos, ope, arch, ΔT)
    in_domain!(phytos, arch, g)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(phytos, ope, arch::Architecture, g::Grids, vel₀, vel½, vel₁, ΔT)
    update_adv_ope!(phytos, ope, arch)

    find_inds!(ope, arch, g, 3, 0)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel₀.u, vel₀.v, vel₀.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 0)

    calc_1st_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    find_inds!(ope, arch, g, 3, 15)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 1)

    calc_2nd_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    find_inds!(ope, arch, g, 3, 15)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel½.u, vel½.v, vel½.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 2)

    calc_3rd_intermediate_coord!(ope, arch, ΔT)
    in_domain!(ope, arch, g, 15)

    find_inds!(ope, arch, g, 3, 15)
    find_vel!(ope, Int.(ope[:,4:6]), arch, g, vel₁.u, vel₁.v, vel₁.w)
    find_xᵈ!(ope, Int.(ope[:,4:6]), arch, g)
    vel_interpolation!(ope, arch, g, 3)

    calc_vel_rk4!(ope, arch)
    calc_coord!(phytos, ope, arch, ΔT)
    in_domain!(phytos, arch, g)
end
