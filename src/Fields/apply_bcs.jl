#####
##### apply flux boundary conditions to tracer tendency `Gc`
##### 

# a positive west flux is associated with an *increase* in `Gc` near the west boundary.
# same for south and bottom.
@kernel function apply_west_bcs_kernel!(Gc, grid, west_bc, iter, ΔT)
    j, k = @index(Global, NTuple)
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[1+grid.Hx, jj, kk] += getbc(west_bc, j, k, iter) * ΔT * Ax(1+grid.Hx, jj, kk, grid) / volume(1+grid.Hx, jj, kk, grid)
end

@kernel function apply_south_bcs_kernel!(Gc, grid, south_bc, iter, ΔT)
    i, k = @index(Global, NTuple)
    ii = i + grid.Hx
    kk = k + grid.Hz
    @inbounds Gc[ii, 1+grid.Hy, kk] += getbc(south_bc, i, k, iter) * ΔT * Ay(ii, 1+grid.Hy, kk, grid) / volume(ii, 1+grid.Hy, kk, grid)
end

@kernel function apply_bottom_bcs_kernel!(Gc, grid, bottom_bc, iter, ΔT)
    i, j = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hy
    @inbounds Gc[ii, jj, grid.Nz+grid.Hz] += getbc(bottom_bc, i, j, iter) * ΔT * Az(ii, jj, grid.Nz+grid.Hz, grid) / volume(ii, jj, grid.Nz+grid.Hz, grid)
end

# a positive east flux is associated with an *decrease* in `Gc` near the east boundary.
# same for north and top.
@kernel function apply_east_bcs_kernel!(Gc, grid, east_bc, iter, ΔT)
    j, k = @index(Global, NTuple)
    jj = j + grid.Hy
    kk = k + grid.Hz
    @inbounds Gc[grid.Nx+grid.Hx, jj, kk] -= getbc(east_bc, j, k, iter) * ΔT * Ax(grid.Nx+grid.Hx, jj, kk, grid) / volume(grid.Nx+grid.Hx, jj, kk, grid)
end

@kernel function apply_north_bcs_kernel!(Gc, grid, north_bc, iter, ΔT)
    i, k = @index(Global, NTuple)
    ii = i + grid.Hx
    kk = k + grid.Hz
    @inbounds Gc[ii, grid.Ny+grid.Hy, kk] -= getbc(north_bc, i, k, iter) * ΔT * Ay(ii, grid.Ny+grid.Hy, kk, grid) / volume(ii, grid.Ny+grid.Hy, kk, grid)
end

@kernel function apply_top_bcs_kernel!(Gc, grid, top_bc, iter, ΔT)
    i, j = @index(Global, NTuple)
    ii = i + grid.Hx
    jj = j + grid.Hy
    @inbounds Gc[ii, jj, 1+grid.Hz] -= getbc(top_bc, i, j, iter) * ΔT * Az(ii, jj, 1+grid.Hz, grid) / volume(ii, jj, 1+grid.Hz, grid)
end

##### apply flux boundary conditions in x direction (west and east)
function apply_west_bcs!(Gc, grid, west_bc, iter, ΔT, arch, dep)
    kernel! = apply_west_bcs_kernel!(device(arch), (16,16), (grid.Ny, grid.Nz))
    event = kernel!(Gc, grid, west_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end
function apply_east_bcs!(Gc, grid, east_bc, iter, ΔT, arch, dep)
    kernel! = apply_east_bcs_kernel!(device(arch), (16,16), (grid.Ny, grid.Nz))
    event = kernel!(Gc, grid, east_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end

##### apply flux boundary conditions in y direction (south and north)
function apply_south_bcs!(Gc, grid, south_bc, iter, ΔT, arch, dep)
    kernel! = apply_south_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Nz))
    event = kernel!(Gc, grid, south_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end
function apply_north_bcs!(Gc, grid, north_bc, iter, ΔT, arch, dep)
    kernel! = apply_north_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Nz))
    event = kernel!(Gc, grid, north_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end

##### apply flux boundary conditions in z direction (bottom and top)
function apply_bottom_bcs!(Gc, grid, bottom_bc, iter, ΔT, arch, dep)
    kernel! = apply_bottom_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny))
    event = kernel!(Gc, grid, bottom_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end
function apply_top_bcs!(Gc, grid, top_bc, iter, ΔT, arch, dep)
    kernel! = apply_top_bcs_kernel!(device(arch), (16,16), (grid.Nx, grid.Ny))
    event = kernel!(Gc, grid, top_bc, iter, ΔT, dependencies=dep)
    wait(device(arch), event)
    return nothing
end


##### apply flux boundary conditions in all three directions, six faces.
function apply_bcs!(Gcs, nuts, grid, iter, ΔT, arch)
    barrier = Event(device(arch))
    events = []
    for name in nut_names
        west_event   =   apply_west_bcs!(Gcs[name].data, grid, nuts[name].bc.west,   iter, ΔT, arch, barrier)
        east_event   =   apply_east_bcs!(Gcs[name].data, grid, nuts[name].bc.east,   iter, ΔT, arch, barrier)
        south_event  =  apply_south_bcs!(Gcs[name].data, grid, nuts[name].bc.south,  iter, ΔT, arch, barrier)
        north_event  =  apply_north_bcs!(Gcs[name].data, grid, nuts[name].bc.north,  iter, ΔT, arch, barrier)
        bottom_event = apply_bottom_bcs!(Gcs[name].data, grid, nuts[name].bc.bottom, iter, ΔT, arch, barrier)
        top_event    =    apply_top_bcs!(Gcs[name].data, grid, nuts[name].bc.top,    iter, ΔT, arch, barrier)
        push!(events, west_event, east_event, south_event, north_event, bottom_event, top_event)
    end

    # filter out ::Nothing
    events = filter(e -> typeof(e) <: Event, events)

    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

# Avoid some computation / memory accesses when no flux boundary conditions will be applied
  apply_west_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing
  apply_east_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing
 apply_north_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing
 apply_south_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing
apply_bottom_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing
   apply_top_bcs!(Gc, grid, ::Nothing, iter, ΔT, arch, dep) = nothing