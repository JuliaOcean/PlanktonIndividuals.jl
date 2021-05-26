##### multi dimensional advection
##### For incompressible model only
##### calculate tendencies in x direction
@kernel function calc_Gcˣ_kernel!(Gc, c, u, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_x(ii, jj, kk, g, u, c, ΔT)
end
function calc_Gcsˣ!(Gcs, nut, u, g::AbstractGrid, ΔT, arch::Architecture)
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
@kernel function calc_Gcʸ_kernel!(Gc, c, v, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_y(ii, jj, kk, g, v, c, ΔT)
end
function calc_Gcsʸ!(Gcs, nut, v, g::AbstractGrid, ΔT, arch::Architecture)
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
@kernel function calc_Gcᶻ_kernel!(Gc, c, w, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds Gc[ii, jj, kk] = adv_flux_z(ii, jj, kk, g, w, c, ΔT)
end
function calc_Gcsᶻ!(Gcs, nut, w, g::AbstractGrid, ΔT, arch::Architecture)
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
@kernel function multi_dim_x_kernel!(ctemp, Gc, c, u, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δx⁺(ii, jj, kk, Gc) - c[ii, jj, kk] * δx⁺(ii, jj, kk, g, Trans_x, u))
end
function multi_dim_x!(nut_temp, Gcs, nut, u, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_x_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut_temp[name].data, Gcs[name].data, nut[name].data, u, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### apply the tendency in y direction to tracer c
@kernel function multi_dim_y_kernel!(ctemp, Gc, c, v, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δy⁺(ii, jj, kk, Gc) - c[ii, jj, kk] * δy⁺(ii, jj, kk, g, Trans_y, v ))
end
function multi_dim_y!(nut_temp, Gcs, nut, v, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_y_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut_temp[name].data, Gcs[name].data, nut[name].data, v, g, ΔT, dependencies=barrier)
        push!(events,event)
    end
    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

##### apply the tendency in z direction to tracer c
@kernel function multi_dim_z_kernel!(ctemp, Gc, c, w, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δz⁺(ii, jj, kk, Gc) - c[ii, jj, kk] * δz⁺(ii, jj, kk, g, Trans_z, w ))
end
function multi_dim_z!(nut_temp, Gcs, nut, w, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_z_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        event = kernel!(nut_temp[name].data, Gcs[name].data, nut[name].data, w, g, ΔT, dependencies=barrier)
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

function nut_advection!(nut, nut_temp, Gcs, vel, g::AbstractGrid, ΔT, arch::Architecture)
    for name in nut_names
        @inbounds nut_temp[name].data .= nut[name].data
    end
    ##### x direction
    calc_Gcsˣ!(Gcs, nut_temp, vel.u.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_x!(nut_temp, Gcs, nut, vel.u.data, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    ##### y direction
    calc_Gcsʸ!(Gcs, nut_temp, vel.v.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_y!(nut_temp, Gcs, nut, vel.v.data, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    ##### z direction
    calc_Gcsᶻ!(Gcs, nut_temp, vel.w.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_z!(nut_temp, Gcs, nut, vel.w.data, g, ΔT, arch)
    fill_halo_nut!(nut_temp, g)

    calc_nut_tendency!(Gcs, nut_temp, nut, ΔT)

end