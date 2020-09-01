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

    @views @. c[:, :, 1:grid.Hz] = c[:, :, grid.Hx+1:grid.Hx+1]       # bottom

    @views @. c[grid.Nx+grid.Hx+1:grid.Nx+2*grid.Hx, :, :] = c[1+grid.Hx:2*grid.Hx, :, :] # east
    @views @. c[:, grid.Ny+grid.Hy+1:grid.Ny+2*grid.Hy, :] = c[:, 1+grid.Hy:2*grid.Hy, :] # north

    @views @. c[:, :, grid.Nz+grid.Hz+1:grid.Nz+2*grid.Hz+1] = c[:, :, grid.Nz+grid.Hz:grid.Nz+grid.Hz] # top

    return nothing
end

@inline function fill_halo!(nut::nutrient_fields, grid)
    fill_halo!(nut.DIC, grid)
    fill_halo!(nut.DOC, grid)
    fill_halo!(nut.POC, grid)
    fill_halo!(nut.NO3, grid)
    fill_halo!(nut.NH4, grid)
    fill_halo!(nut.DON, grid)
    fill_halo!(nut.PON, grid)
    fill_halo!(nut.PO4, grid)
    fill_halo!(nut.DOP, grid)
    fill_halo!(nut.POP, grid)
end

@inline function zero_halo!(nut::nutrient_fields, grid)
    zero_halo!(nut.DIC, grid)
    zero_halo!(nut.DOC, grid)
    zero_halo!(nut.POC, grid)
    zero_halo!(nut.NO3, grid)
    zero_halo!(nut.NH4, grid)
    zero_halo!(nut.DON, grid)
    zero_halo!(nut.PON, grid)
    zero_halo!(nut.PO4, grid)
    zero_halo!(nut.DOP, grid)
    zero_halo!(nut.POP, grid)
end
