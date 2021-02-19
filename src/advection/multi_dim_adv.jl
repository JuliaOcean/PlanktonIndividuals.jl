##### multi dimensional advection
##### For incompressible model only
##### calculate tendencies in x direction
@kernel function calc_Gcˣ_kernel!(Gc, c, u, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_x(ii, jj, kk, g, u, c, ΔT) / g.V
end
function calc_Gcsˣ!(Gcs, nut, u, g::Grids, ΔT, arch::Architecture)
    kernel! = calc_Gcˣ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(Gcs[name].data, nut[name].data, u, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### calculate tendencies in y direction
@kernel function calc_Gcʸ_kernel!(Gc, c, v, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_y(ii, jj, kk, g, v, c, ΔT) / g.V
end
function calc_Gcsʸ!(Gcs, nut, v, g::Grids, ΔT, arch::Architecture)
    kernel! = calc_Gcʸ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(Gcs[name].data, nut[name].data, v, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### calculate tendencies in z direction
@kernel function calc_Gcᶻ_kernel!(Gc, c, w, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_z(ii, jj, kk, g, w, c, ΔT) / g.V
end
function calc_Gcsᶻ!(Gcs, nut, w, g::Grids, ΔT, arch::Architecture)
    kernel! = calc_Gcᶻ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(Gcs[name].data, nut[name].data, w, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### apply the tendency in x direction to tracer c
@kernel function multi_dim_x_kernel!(ctemp, Gc, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT * δx⁺(ii, jj, kk, Gc)
end
function multi_dim_x!(nut, Gcs, g::Grids, ΔT, arch::Architecture)
    kernel! = multi_dim_x_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut[name].data, Gcs[name].data, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### apply the tendency in y direction to tracer c
@kernel function multi_dim_y_kernel!(ctemp, Gc, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT * δy⁺(ii, jj, kk, Gc)
end
function multi_dim_y!(nut, Gcs, g::Grids, ΔT, arch::Architecture)
    kernel! = multi_dim_y_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut[name].data, Gcs[name].data, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### apply the tendency in z direction to tracer c
@kernel function multi_dim_z_kernel!(ctemp, Gc, g::Grids, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT * δz⁺(ii, jj, kk, Gc)
end
function multi_dim_z!(nut, Gcs, g::Grids, ΔT, arch::Architecture)
    kernel! = multi_dim_z_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut[name].data, Gcs[name].data, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

function calc_nut_tendency!(a, b, c, ΔT)
    for name in nut_names
        @inbounds a[name].data .= (b[name].data .- c[name].data) ./ ΔT
    end
end

function nut_advection!(nut, nut_temp, Gcs, vel, g::Grids, ΔT, arch::Architecture)
    for name in nut_names
        @inbounds nut_temp[name].data .= nut[name].data
    end
    ##### x direction
    calc_Gcsˣ!(Gcs, nut_temp, vel.u.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_x!(nut_temp, Gcs, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    ##### y direction
    calc_Gcsʸ!(Gcs, nut_temp, vel.v.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_y!(nut_temp, Gcs, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    ##### z direction
    calc_Gcsᶻ!(Gcs, nut_temp, vel.w.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_z!(nut_temp, Gcs, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    calc_nut_tendency!(Gcs, nut_temp, nut, ΔT)

end