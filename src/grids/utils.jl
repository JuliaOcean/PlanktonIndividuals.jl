#####
##### replace the storage place of the grid information based on the architecture
#####
function replace_grid_storage(arch::Architecture, grid::RegularLatLonGrid{TX, TY, TZ}) where {TX, TY, TZ}
    dxC = grid.dxC |> array_type(arch)
    dyC = grid.dyC |> array_type(arch)
    dxF = grid.dxF |> array_type(arch)
    dyF = grid.dyF |> array_type(arch)
    Ax  = grid.Ax  |> array_type(arch)
    Ay  = grid.Ay  |> array_type(arch)
    Az  = grid.Az  |> array_type(arch)
    Vol = grid.Vol |> array_type(arch)

    return RegularLatLonGrid{TX, TY, TZ, typeof(grid.xF)}(
        grid.xC, grid.yC, grid.zC, grid.xF, grid.yF, grid.zF, grid.Δx, grid.Δy, grid.Δz,
        dxC, dyC, dxF, dyF, Ax, Ay, Az, Vol, grid.Nx, grid.Ny, grid.Nz, grid.Hx, grid.Hy, grid.Hz)
end

function replace_grid_storage(arch::Architecture, grid::VerticallyStretchedLatLonGrid{TX, TY, TZ}) where {TX, TY, TZ}
    dxC = grid.dxC |> array_type(arch)
    dyC = grid.dyC |> array_type(arch)
    dzC = grid.dzC |> array_type(arch)
    dxF = grid.dxF |> array_type(arch)
    dyF = grid.dyF |> array_type(arch)
    dzF = grid.dzF |> array_type(arch)
    Ax  = grid.Ax  |> array_type(arch)
    Ay  = grid.Ay  |> array_type(arch)
    Az  = grid.Az  |> array_type(arch)
    Vol = grid.Vol |> array_type(arch)

    return VerticallyStretchedLatLonGrid{TX, TY, TZ, typeof(grid.xF)}(
        grid.xC, grid.yC, grid.zC, grid.xF, grid.yF, grid.zF, grid.Δx, grid.Δy, dxC, dyC, dzC,
        dxF, dyF, dzF, Ax, Ay, Az, Vol, grid.Nx, grid.Ny, grid.Nz, grid.Hx, grid.Hy, grid.Hz)
end

function replace_grid_storage(arch::Architecture, grid::RegularRectilinearGrid{TX, TY, TZ}) where {TX, TY, TZ}
    return grid
end

#####
##### adapt the grid struct to GPU
#####

Adapt.adapt_structure(to, grid::RegularLatLonGrid{TX, TY, TZ}) where {TX, TY, TZ} =
    RegularLatLonGrid{TX, TY, TZ, typeof(grid.xF)}(
        grid.xC, grid.yC, grid.zC,
        grid.xF, grid.yF, grid.zF,
        grid.Δx, grid.Δy, grid.Δz,
        Adapt.adapt(to, grid.dxC),
        Adapt.adapt(to, grid.dyC),
        Adapt.adapt(to, grid.dxF),
        Adapt.adapt(to, grid.dyF),
        Adapt.adapt(to, grid.Ax),
        Adapt.adapt(to, grid.Ay),
        Adapt.adapt(to, grid.Az),
        Adapt.adapt(to, grid.Vol),
        grid.Nx, grid.Ny, grid.Nz,
        grid.Hx, grid.Hy, grid.Hz)

Adapt.adapt_structure(to, grid::VerticallyStretchedLatLonGrid{TX, TY, TZ}) where {TX, TY, TZ} =
    VerticallyStretchedLatLonGrid{TX, TY, TZ, typeof(grid.xF)}(
        grid.xC, grid.yC, grid.zC,
        grid.xF, grid.yF, grid.zF,
        grid.Δx, grid.Δy,
        Adapt.adapt(to, grid.dxC),
        Adapt.adapt(to, grid.dyC),
        Adapt.adapt(to, grid.dzC),
        Adapt.adapt(to, grid.dxF),
        Adapt.adapt(to, grid.dyF),
        Adapt.adapt(to, grid.dzF),
        Adapt.adapt(to, grid.Ax),
        Adapt.adapt(to, grid.Ay),
        Adapt.adapt(to, grid.Az),
        Adapt.adapt(to, grid.Vol),
        grid.Nx, grid.Ny, grid.Nz,
        grid.Hx, grid.Hy, grid.Hz)