##### multi dimensional advection
##### For incompressible model only

##### calculate the divergence of the flux of tracer q only advected by u or v or w
@inline div_flux_x(i, j, k, g::grids, u, q, ΔT) = 1/g.V * δx⁺(i, j, k, g::grids, adv_flux_x, u, q, ΔT)
@inline div_flux_y(i, j, k, g::grids, v, q, ΔT) = 1/g.V * δy⁺(i, j, k, g::grids, adv_flux_y, v, q, ΔT)
@inline div_flux_z(i, j, k, g::grids, w, q, ΔT) = 1/g.V * δz⁺(i, j, k, g::grids, adv_flux_z, w, q, ΔT)

##### apply the tendency of by multi-dimensional advection for tracer q
@kernel function calc_qˣ!(qtemp, g::grids, u, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_x(ii, jj, kk, g::grids, u, q, ΔT)
end
@kernel function calc_qʸ!(qtemp, g::grids, v, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_y(ii, jj, kk, g::grids, v, q, ΔT)
end
@kernel function calc_qᶻ!(qtemp, g::grids, w, q, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds qtemp[ii, jj, kk] = q[ii, jj, kk] - ΔT * div_flux_z(ii, jj, kk, g::grids, w, q, ΔT)
end

function nut_advection!(nutₜ, arch::Architecture, g, nutrients, vel, ΔT)
    calc_qˣ_kernel! = calc_qˣ(device(arch), (g.Nx, g.Ny, g.Nz))
    calc_qʸ_kernel! = calc_qʸ(device(arch), (g.Nx, g.Ny, g.Nz))
    calc_qᶻ_kernel! = calc_qᶻ(device(arch), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))

    DIC_event_x = calc_qˣ_kernel!(nutₜ.DIC, g, vel.u, nutrients.DIC, ΔT, dependencies=barrier)
    NH4_event_x = calc_qˣ_kernel!(nutₜ.NH4, g, vel.u, nutrients.NH4, ΔT, dependencies=barrier)
    NO3_event_x = calc_qˣ_kernel!(nutₜ.NO3, g, vel.u, nutrients.NO3, ΔT, dependencies=barrier)
    PO4_event_x = calc_qˣ_kernel!(nutₜ.PO4, g, vel.u, nutrients.PO4, ΔT, dependencies=barrier)
    DOC_event_x = calc_qˣ_kernel!(nutₜ.DOC, g, vel.u, nutrients.DOC, ΔT, dependencies=barrier)
    POC_event_x = calc_qˣ_kernel!(nutₜ.POC, g, vel.u, nutrients.POC, ΔT, dependencies=barrier)
    DON_event_x = calc_qˣ_kernel!(nutₜ.DON, g, vel.u, nutrients.DON, ΔT, dependencies=barrier)
    PON_event_x = calc_qˣ_kernel!(nutₜ.PON, g, vel.u, nutrients.PON, ΔT, dependencies=barrier)
    DOP_event_x = calc_qˣ_kernel!(nutₜ.DOP, g, vel.u, nutrients.DOP, ΔT, dependencies=barrier)
    POP_event_x = calc_qˣ_kernel!(nutₜ.POP, g, vel.u, nutrients.POP, ΔT, dependencies=barrier)

    events_x = [DIC_events_x, NH4_events_x, NO3_events_x, PO4_events_x,
                DOC_events_x, POC_events_x, DON_events_x, PON_events_x,
                DOP_events_x, POP_events_x]

    wait(device(arch), MultiEvent(Tuple(events_x)))

    DIC_event_y = calc_ʸq_kernel!(nutₜ.DIC, g, vel.v, nutₜ.DIC, ΔT, dependencies=barrier)
    NH4_event_y = calc_ʸq_kernel!(nutₜ.NH4, g, vel.v, nutₜ.NH4, ΔT, dependencies=barrier)
    NO3_event_y = calc_ʸq_kernel!(nutₜ.NO3, g, vel.v, nutₜ.NO3, ΔT, dependencies=barrier)
    PO4_event_y = calc_ʸq_kernel!(nutₜ.PO4, g, vel.v, nutₜ.PO4, ΔT, dependencies=barrier)
    DOC_event_y = calc_ʸq_kernel!(nutₜ.DOC, g, vel.v, nutₜ.DOC, ΔT, dependencies=barrier)
    POC_event_y = calc_ʸq_kernel!(nutₜ.POC, g, vel.v, nutₜ.POC, ΔT, dependencies=barrier)
    DON_event_y = calc_ʸq_kernel!(nutₜ.DON, g, vel.v, nutₜ.DON, ΔT, dependencies=barrier)
    PON_event_y = calc_ʸq_kernel!(nutₜ.PON, g, vel.v, nutₜ.PON, ΔT, dependencies=barrier)
    DOP_event_y = calc_ʸq_kernel!(nutₜ.DOP, g, vel.v, nutₜ.DOP, ΔT, dependencies=barrier)
    POP_event_y = calc_qʸ_kernel!(nutₜ.POP, g, vel.v, nutₜ.POP, ΔT, dependencies=barrier)

    events_y = [DIC_events_y, NH4_events_y, NO3_events_y, PO4_events_y,
                DOC_events_y, POC_events_y, DON_events_y, PON_events_y,
                DOP_events_y, POP_events_y]

    wait(device(arch), MultiEvent(Tuple(events_y)))

    DIC_event_z = calc_qᶻ_kernel!(nutₜ.DIC, g, vel.w, nutₜ.DIC, ΔT, dependencies=barrier)
    NH4_event_z = calc_qᶻ_kernel!(nutₜ.NH4, g, vel.w, nutₜ.NH4, ΔT, dependencies=barrier)
    NO3_event_z = calc_qᶻ_kernel!(nutₜ.NO3, g, vel.w, nutₜ.NO3, ΔT, dependencies=barrier)
    PO4_event_z = calc_qᶻ_kernel!(nutₜ.PO4, g, vel.w, nutₜ.PO4, ΔT, dependencies=barrier)
    DOC_event_z = calc_qᶻ_kernel!(nutₜ.DOC, g, vel.w, nutₜ.DOC, ΔT, dependencies=barrier)
    POC_event_z = calc_qᶻ_kernel!(nutₜ.POC, g, vel.w, nutₜ.POC, ΔT, dependencies=barrier)
    DON_event_z = calc_qᶻ_kernel!(nutₜ.DON, g, vel.w, nutₜ.DON, ΔT, dependencies=barrier)
    PON_event_z = calc_qᶻ_kernel!(nutₜ.PON, g, vel.w, nutₜ.PON, ΔT, dependencies=barrier)
    DOP_event_z = calc_qᶻ_kernel!(nutₜ.DOP, g, vel.w, nutₜ.DOP, ΔT, dependencies=barrier)
    POP_event_z = calc_qᶻ_kernel!(nutₜ.POP, g, vel.w, nutₜ.POP, ΔT, dependencies=barrier)

    events_z = [DIC_events_z, NH4_events_z, NO3_events_z, PO4_events_z,
                DOC_events_z, POC_events_z, DON_events_z, PON_events_z,
                DOP_events_z, POP_events_z]

    wait(device(arch), MultiEvent(Tuple(events_z)))

    return nothing
end

