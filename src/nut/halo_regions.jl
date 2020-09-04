@inline function zero_halo!(c::AbstractArray, grid)
    @views @. c[1:grid.Hx, :, :] = 0 # west
    @views @. c[:, 1:grid.Hy, :] = 0 # south
    @views @. c[:, :, 1:grid.Hz] = 0 # bottom
    @views @. c[grid.Nx+grid.Hx+1:grid.Nx+2*grid.Hx, :, :] = 0   # east
    @views @. c[:, grid.Ny+grid.Hy+1:grid.Ny+2*grid.Hy, :] = 0   # north
    @views @. c[:, :, grid.Nz+grid.Hz+1:grid.Nz+2*grid.Hz+1] = 0 # top + free surface
    return nothing
end

@inline function fill_halo!(c::AbstractArray, grid)
    @views @. c[1:grid.Hx, :, :] = c[grid.Nx+1:grid.Nx+grid.Hx, :, :] # west
    @views @. c[:, 1:grid.Hy, :] = c[:, grid.Ny+1:grid.Ny+grid.Hy, :] # south

    @views @. c[:, :, 1:grid.Hz] = c[:, :, grid.Hz+1:grid.Hz+1]       # bottom

    @views @. c[grid.Nx+grid.Hx+1:grid.Nx+2*grid.Hx, :, :] = c[1+grid.Hx:2*grid.Hx, :, :] # east
    @views @. c[:, grid.Ny+grid.Hy+1:grid.Ny+2*grid.Hy, :] = c[:, 1+grid.Hy:2*grid.Hy, :] # north

    @views @. c[:, :, grid.Nz+grid.Hz+1:grid.Nz+2*grid.Hz+1] = c[:, :, grid.Nz+grid.Hz:grid.Nz+grid.Hz] # top

    return nothing
end

@inline function fill_halo!(nuts::NamedTuple, grid)
    for nut in nuts
        fill_halo!(nut.data, grid)
    end

    return nothing
end

@inline function zero_halo!(nuts::NamedTuple, grid)
    for nut in nuts
        zero_halo!(nut.data, grid)
    end

    return nothing
end
