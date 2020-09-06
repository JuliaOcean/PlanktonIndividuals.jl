##### update coordinates of each individual using Explicit Euler (aka Euler Forward) integration
function plankton_advection!(phytos, arch::Architecture, g::Grids, vel, ΔT)
    op_array = adv_op_array_setup(phytos, arch)
    find_inds!(op_array, arch, g)
    vel₁_interpolation!(op_array, arch, g, vel.u, vel.v, vel.w)
    calc_coord!(op_array, phytos, arch, g, ΔT)
end

##### update coordinates of each individual using RK4 integration
function plankton_advectionRK4!(phytos, arch::Architecture, g::Grids, vel₀, vel₁, ΔT)
    vel½ = (u = (vel₀.u .+ vel₁.u) .* 0.5, v = (vel₀.v .+ vel₁.v) .* 0.5, w = (vel₀.w .+ vel₁.w) .* 0.5)
    op_array = adv_op_array_setup(phytos, arch)

    find_inds!(op_array, arch, g)
    vel₁_interpolation!(op_array, arch, g, vel₀.u, vel₀.v, vel₀.w)
    calc_intermediate1_coord!(op_array, phytos, arch, g, ΔT)

    update_inds!(op_array, arch, g)
    vel₂_interpolation!(op_array, arch, g, vel½.u, vel½.v, vel½.w)
    calc_intermediate2_coord!(op_array, phytos, arch, g, ΔT)

    update_inds!(op_array, arch, g)
    vel₃_interpolation!(op_array, arch, g, vel½.u, vel½.v, vel½.w)
    calc_intermediate3_coord!(op_array, phytos, arch, g, ΔT)

    update_inds!(op_array, arch, g)
    vel₄_interpolation!(op_array, arch, g, vel₁.u, vel₁.v, vel₁.w)
    calc_vel_rk4!(op_array, arch)
    calc_coord!(op_array, phytos, arch, g, ΔT)
end
