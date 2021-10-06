#####
##### apply flux boundary conditions to tracer tendency `Gc`
##### 

# a positive west flux is associated with an *increase* in `Gc` near the west boundary.
# same for south and bottom.
@inline function apply_west_bc!(Gc, west_flux, j, k, iter, grid, ΔT)
    @inbounds Gc[1+grid.Hx, j, k] += getbc(west_flux, j, k, iter) * ΔT * Ax(1+grid.Hx, j, k, grid) / volume(1+grid.Hx, j, k, grid)
    return nothing
end

@inline function apply_south_bc!(Gc, south_flux, i, k, iter, grid, ΔT)
    @inbounds Gc[i, 1+grid.Hy, k] += getbc(south_flux, i, k, iter) * ΔT * Ay(i, 1+grid.Hy, k, grid) / volume(i, 1+grid.Hy, k, grid)
    return nothing
end

@inline function apply_bottom_bc!(Gc, bottom_flux, i, j, iter, grid, ΔT)
    @inbounds Gc[i, j, 1+grid.Hz] += getbc(bottom_flux, i, j, iter) * ΔT * Az(i, j, 1+grid.Hz, grid) / volume(i, j, 1+grid.Hz, grid)
    return nothing
end

# a positive east flux is associated with an *decrease* in `Gc` near the east boundary.
# same for north and top.
@inline function apply_east_bc!(Gc, east_flux, j, k, iter, grid, ΔT)
    @inbounds Gc[grid.Nx+grid.Hx, j, k] += getbc(east_flux, j, k, iter) * ΔT * Ax(grid.Nx+grid.Hx, j, k, grid) / volume(grid.Nx+grid.Hx, j, k, grid)
    return nothing
end

@inline function apply_north_bc!(Gc, north_flux, i, k, iter, grid, ΔT)
    @inbounds Gc[i, grid.Ny+grid.Hy, k] += getbc(north_flux, i, k, iter) * ΔT * Ay(i, grid.Ny+grid.Hy, k, grid) / volume(i, grid.Ny+grid.Hy, k, grid)
    return nothing
end

@inline function apply_top_bc!(Gc, top_flux, i, j, iter, grid, ΔT)
    @inbounds Gc[i, j, grid.Nz+grid.Hz] += getbc(top_flux, i, j, iter) * ΔT * Az(i, j, grid.Nz+grid.Hz, grid) / volume(i, j, grid.Nz+grid.Hz, grid)
    return nothing
end

##### apply flux boundary conditions in x direction (west and east)
@kernel function apply_x_bcs_kernel!(Gc, grid, west_bc, east_bc, iter, ΔT)
    j, k = @index(Global, NTuple)
    jj = j + grid.Hy
    kk = k + grid.Hz
    apply_west_bc!(Gc, west_bc, jj, kk, iter, grid, ΔT)
    apply_east_bc!(Gc, east_bc, jj, kk, iter, grid, ΔT)
end
function apply_x_bcs!(Gc, grid, west_bc, east_bc, iter, ΔT, arch, dep)
    kernel! = apply_x_bcs_kernel!(device(arch), (16,16), (grid.Ny, grid.Nz))
    event = kernel!(Gc, grid, west_bc, east_bc, iter, ΔT, dependencies=dep )
    wait(device(arch), event)
    return nothing
end

##### apply flux boundary conditions in y direction (south and north)
@kernel function apply_y_bcs_kernel!(Gc, grid, south_bc, north_bc, iter, ΔT)
    i, k = @index(Global, NTuple)
    ii = i + grid.Hx
    kk = k + grid.Hz
    apply_south_bc!(Gc, south_bc, ii, kk, iter, grid, ΔT)
    apply_north_bc!(Gc, north_bc, ii, kk, iter, grid, ΔT)
end
function apply_y_bcs!(Gc, grid, south_bc, north_bc, iter, ΔT, arch, dep)
    kernel! = apply_y_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Nz))
    event = kernel!(Gc, grid, south_bc, north_bc, iter, ΔT, dependencies=dep )
    wait(device(arch), event)
    return nothing
end

##### apply flux boundary conditions in z direction (bottom and top)
@kernel function apply_z_bcs_kernel!(Gc, grid, bottom_bc, top_bc, iter, ΔT)
    i, j = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hx
    apply_bottom_bc!(Gc, bottom_bc, ii, jj, iter, grid, ΔT)
    apply_top_bc!(Gc, top_bc, ii, jj, iter, grid, ΔT)
end
function apply_z_bcs!(Gc, grid, bottom_bc, top_bc, iter, ΔT, arch, dep)
    kernel! = apply_z_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny))
    event = kernel!(Gc, grid, bottom_bc, top_bc, iter, ΔT, dependencies=dep )
    wait(device(arch), event)
    return nothing
end


##### apply flux boundary conditions in all three directions, six faces.
function apply_bcs!(Gcs, nuts, grid, iter, ΔT, arch)
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        x_event = apply_x_bcs!(Gcs[name].data, grid, nuts[name].bc.x.left, nuts[name].bc.x.right, iter, ΔT, arch, barrier)
        y_event = apply_y_bcs!(Gcs[name].data, grid, nuts[name].bc.y.left, nuts[name].bc.y.right, iter, ΔT, arch, barrier)
        z_event = apply_z_bcs!(Gcs[name].data, grid, nuts[name].bc.z.left, nuts[name].bc.z.right, iter, ΔT, arch, barrier)
        push!(events, x_event, y_event, z_event)
    end

    # filter out ::Nothing
    events = filter(e -> typeof(e) <: Event, events)

    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

# Avoid some computation / memory accesses when no flux boundary conditions will be applied
apply_x_bcs!(Gc, grid, ::Nothing, ::Nothing, iter, ΔT, arch, dep) = nothing
apply_y_bcs!(Gc, grid, ::Nothing, ::Nothing, iter, ΔT, arch, dep) = nothing
apply_z_bcs!(Gc, grid, ::Nothing, ::Nothing, iter, ΔT, arch, dep) = nothing
