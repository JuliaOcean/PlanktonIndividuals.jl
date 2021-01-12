@inline function zero_halo_west!(c::AbstractArray, grid)
    @views @. c[1:grid.Hx, :, :] = 0 # west
    return nothing
end
@inline function zero_halo_east!(c::AbstractArray, grid)
    @views @. c[grid.Nx+grid.Hx+1:grid.Nx+2*grid.Hx, :, :] = 0   # east
    return nothing
end
@inline function zero_halo_south!(c::AbstractArray, grid)
    @views @. c[:, 1:grid.Hy, :] = 0 # south
    return nothing
end
@inline function zero_halo_north!(c::AbstractArray, grid)
    @views @. c[:, grid.Ny+grid.Hy+1:grid.Ny+2*grid.Hy, :] = 0   # north
    return nothing
end
@inline function zero_halo_bottom!(c::AbstractArray, grid)
    @views @. c[:, :, 1:grid.Hz] = 0 # bottom
    return nothing
end
@inline function zero_halo_top!(c::AbstractArray, grid)
    @views @. c[:, :, grid.Nz+grid.Hz+1:grid.Nz+2*grid.Hz] = 0   # top
    return nothing
end

@inline function fill_halo_west!(c::AbstractArray, grid)
    @views @. c[1:grid.Hx, :, :] = c[grid.Nx+1:grid.Nx+grid.Hx, :, :] # west
    return nothing
end
@inline function fill_halo_east!(c::AbstractArray, grid)
    @views @. c[grid.Nx+grid.Hx+1:grid.Nx+2*grid.Hx, :, :] = c[1+grid.Hx:2*grid.Hx, :, :] # east
    return nothing
end
@inline function fill_halo_south!(c::AbstractArray, grid)
    @views @. c[:, 1:grid.Hy, :] = c[:, grid.Ny+1:grid.Ny+grid.Hy, :] # south
    return nothing
end
@inline function fill_halo_north!(c::AbstractArray, grid)
    @views @. c[:, grid.Ny+grid.Hy+1:grid.Ny+2*grid.Hy, :] = c[:, 1+grid.Hy:2*grid.Hy, :] # north
    return nothing
end
@inline function fill_halo_bottom!(c::AbstractArray, grid)
    @views @. c[:, :, 1:grid.Hz] = c[:, :, grid.Hz+1:grid.Hz+1] # bottom
    return nothing
end
@inline function fill_halo_top!(c::AbstractArray, grid)
    @views @. c[:, :, grid.Nz+grid.Hz+1:grid.Nz+2*grid.Hz] = c[:, :, grid.Nz+grid.Hz:grid.Nz+grid.Hz] # top
    return nothing
end

@inline function fill_halo_nut!(nuts::NamedTuple, grid)
    for nut in nuts
        fill_halo_west!(nut.data, grid)
        fill_halo_east!(nut.data, grid)
        fill_halo_south!(nut.data, grid)
        fill_halo_north!(nut.data, grid)
        fill_halo_bottom!(nut.data, grid)
        fill_halo_top!(nut.data, grid)
    end
    return nothing
end
@inline function fill_halo_vel!(vels::NamedTuple, grid)
    for vel in vels
        fill_halo_west!(vel.data, grid)
        fill_halo_east!(vel.data, grid)
        fill_halo_south!(vel.data, grid)
        fill_halo_north!(vel.data, grid)
        zero_halo_bottom!(vel.data, grid)
        zero_halo_top!(vel.data, grid)
    end
    return nothing
end

@inline function zero_halo!(nuts::NamedTuple, grid)
    for nut in nuts
        zero_halo!(nut.data, grid)
    end

    return nothing
end
