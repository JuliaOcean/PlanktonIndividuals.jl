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
function calc_Gcsˣ!(Gcs, tracers, u, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = calc_Gcˣ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(Gcs[name].data, tracers[name].data, u, g, ΔT)
    end
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
function calc_Gcsʸ!(Gcs, tracers, v, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = calc_Gcʸ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(Gcs[name].data, tracers[name].data, v, g, ΔT)
    end
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
function calc_Gcsᶻ!(Gcs, tracers, w, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = calc_Gcᶻ_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(Gcs[name].data, tracers[name].data, w, g, ΔT)
    end
    return nothing
end

##### apply the tendency in x direction to tracer c
@kernel function multi_dim_x_kernel!(ctemp, Gc, c, u, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δx⁺(ii, jj, kk, Gc, g) - c[ii, jj, kk] * δx⁺(ii, jj, kk, g, Trans_x, u))
end
function multi_dim_x!(tracer_temp, Gcs, tracers, u, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_x_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(tracer_temp[name].data, Gcs[name].data, tracers[name].data, u, g, ΔT)
    end
    return nothing
end

##### apply the tendency in y direction to tracer c
@kernel function multi_dim_y_kernel!(ctemp, Gc, c, v, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δy⁺(ii, jj, kk, Gc, g) - c[ii, jj, kk] * δy⁺(ii, jj, kk, g, Trans_y, v ))
end
function multi_dim_y!(tracer_temp, Gcs, tracers, v, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_y_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(tracer_temp[name].data, Gcs[name].data, tracers[name].data, v, g, ΔT)
    end
    return nothing
end

##### apply the tendency in z direction to tracer c
@kernel function multi_dim_z_kernel!(ctemp, Gc, c, w, g::AbstractGrid, ΔT)
    i, j, k = @index(Global, NTuple)
    ### offset index for halo points
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz
    @inbounds ctemp[ii, jj, kk] -= ΔT / volume(ii, jj, kk, g) * (δz⁺(ii, jj, kk, Gc, g) - c[ii, jj, kk] * δz⁺(ii, jj, kk, g, Trans_z, w ))
end
function multi_dim_z!(tracer_temp, Gcs, tracers, w, g::AbstractGrid, ΔT, arch::Architecture)
    kernel! = multi_dim_z_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    for name in tracer_names
        kernel!(tracer_temp[name].data, Gcs[name].data, tracers[name].data, w, g, ΔT)
    end
    return nothing
end

function calc_tracer_tendency!(a, b, c)
    for name in tracer_names
        @inbounds a[name].data .= b[name].data .- c[name].data
    end
end

function tracer_advection!(tracer, tracer_temp, Gcs, vel, g::AbstractGrid, ΔT, arch::Architecture)
    for name in tracer_names
        @inbounds tracer_temp[name].data .= tracer[name].data
    end
    ##### x direction
    calc_Gcsˣ!(Gcs, tracer_temp, vel.u.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_x!(tracer_temp, Gcs, tracer, vel.u.data, g, ΔT, arch)
    fill_halo_tracer!(tracer_temp, g)

    ##### y direction
    calc_Gcsʸ!(Gcs, tracer_temp, vel.v.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_y!(tracer_temp, Gcs, tracer, vel.v.data, g, ΔT, arch)
    fill_halo_tracer!(tracer_temp, g)

    ##### z direction
    calc_Gcsᶻ!(Gcs, tracer_temp, vel.w.data, g, ΔT, arch)
    fill_halo_Gcs!(Gcs, g)
    multi_dim_z!(tracer_temp, Gcs, tracer, vel.w.data, g, ΔT, arch)
    fill_halo_tracer!(tracer_temp, g)

    calc_tracer_tendency!(Gcs, tracer_temp, tracer)

end
