@inline function fill_halo_horizontally_periodic(c, g::Grids)
    @views @. c[1:g.Hx, :, :] = c[g.Nx+1:g.Nx+g.Hx, :, :] # west
    @views @. c[:, 1:g.Hy, :] = c[:, g.Ny+1:g.Ny+g.Hy, :] # south
    @views @. c[g.Nx+g.Hx+1:g.Nx+2*g.Hx, :, :] = c[1+g.Hx:2*g.Hx, :, :] # east
    @views @. c[:, g.Ny+g.Hy+1:g.Ny+2*g.Hy, :] = c[:, 1+g.Hy:2*g.Hy, :] # north
    return nothing
end
@inline function fill_halo_vertically_bounded(c, g::Grids)
    @views @. c[:, :, 1:g.Hz] = c[:, :, g.Hz+1:g.Hz+1] # bottom
    @views @. c[:, :, g.Nz+g.Hz+1:g.Nz+2*g.Hz] = c[:, :, g.Nz+g.Hz:g.Nz+g.Hz] # top
    return nothing
end

@inline function fill_halo_nut!(nuts::NamedTuple, g::Grids)
    for nut in nuts
        fill_halo_horizontally_periodic(nut.data, g)
        fill_halo_vertically_bounded(nut.data, g)
    end
    return nothing
end
@inline function fill_halo_vel!(vels::NamedTuple, g::Grids)
    for vel in vels
        fill_halo_horizontally_periodic(vel.data, g)
    end
    fill_halo_vertically_bounded(vels.u.data, g)
    fill_halo_vertically_bounded(vels.v.data, g)
    @views @. vels.w.data[:, :, 1:g.Hz] = vels.w.data[:, :, g.Hz+1:g.Hz+1] # bottom
    @views @. vels.w.data[:, :, g.Nz+g.Hz+2:g.Nz+2*g.Hz] = vels.w.data[:, :, g.Nz+g.Hz+1:g.Nz+g.Hz+1] # top
    return nothing
end

@inline function fill_halo_Gcs!(nuts::NamedTuple, g::Grids)
    for nut in nuts
        @views @. nut.data[g.Nx+g.Hx+1:g.Nx+2*g.Hx, :, :] = nut.data[1+g.Hx:2*g.Hx, :, :] # east
        @views @. nut.data[:, g.Ny+g.Hy+1:g.Ny+2*g.Hy, :] = nut.data[:, 1+g.Hy:2*g.Hy, :] # north
        @views @. nut.data[:, :, g.Nz+g.Hz+1:g.Nz+2*g.Hz] = 0.0 # top
    end
    return nothing
end