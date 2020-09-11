##### multi dimensional advection
##### For incompressible model only

##### calculate the divergence of the flux of tracer q only advected by u or v or w
@inline div_flux_x(i, j, k, g::Grids, u, q, ΔT) = 1/g.V * δx⁺(i, j, k, g, adv_flux_x, u, q, ΔT)
@inline div_flux_y(i, j, k, g::Grids, v, q, ΔT) = 1/g.V * δy⁺(i, j, k, g, adv_flux_y, v, q, ΔT)
@inline div_flux_z(i, j, k, g::Grids, w, q, ΔT) = 1/g.V * δz⁺(i, j, k, g, adv_flux_z, w, q, ΔT)

##### apply the tendency of by multi-dimensional advection for tracer q
@kernel function calc_qˣ!(qtemp, g::Grids, u, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_x(ii, jj, kk, g, u, q, ΔT)
end
@kernel function calc_qʸ!(qtemp, g::Grids, v, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_y(ii, jj, kk, g, v, q, ΔT)
end
@kernel function calc_qᶻ!(qtemp, g::Grids, w, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_z(ii, jj, kk, g, w, q, ΔT)
end

function nut_advectionˣ!(nutₜ, arch::Architecture, g, nutrients, u, ΔT)
    calc_qˣ_kernel! = calc_qˣ!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    events_x = []
    for name in nut_names
        event = calc_qˣ_kernel!(nutₜ[name].data, g, u, nutrients[name].data, ΔT, dependencies=barrier)
        push!(events_x,event)
    end

    wait(device(arch), MultiEvent(Tuple(events_x)))

    return nothing
end

function nut_advectionʸ!(nutₜ, arch::Architecture, g, nutrients, v, ΔT)
    calc_qʸ_kernel! = calc_qʸ!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    events_y = []
    for name in nut_names
        event = calc_qʸ_kernel!(nutₜ[name].data, g, v, nutrients[name].data, ΔT, dependencies=barrier)
        push!(events_y,event)
    end

    wait(device(arch), MultiEvent(Tuple(events_y)))

    return nothing
end

function nut_advectionᶻ!(nutₜ, arch::Architecture, g, nutrients, w, ΔT)
    calc_qᶻ_kernel! = calc_qᶻ!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    events_z = []
    for name in nut_names
        event = calc_qᶻ_kernel!(nutₜ[name].data, g, w, nutrients[name].data, ΔT, dependencies=barrier)
        push!(events_z,event)
    end

    wait(device(arch), MultiEvent(Tuple(events_z)))

    return nothing
end

function nut_advection!(Gc, arch::Architecture, g, nut, nut₁, nut₂, nut₃, vel, ΔT)
    nut_advectionˣ!(nut₁, arch::Architecture, g, nut, vel.u, ΔT)

    fill_halo!(nut₁, g)

    nut_advectionʸ!(nut₂, arch::Architecture, g, nut₁, vel.v, ΔT)

    fill_halo!(nut₂, g)

    nut_advectionᶻ!(nut₃, arch::Architecture, g, nut₂, vel.w, ΔT)

    sub_nut_tendency!(Gc, nut₃, nut)

end

