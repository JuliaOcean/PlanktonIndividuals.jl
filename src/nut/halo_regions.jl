@inline function fill_halo_nut!(nuts::NamedTuple, g::Grids)
    for nut in nuts
        @views @. nut.data[1:g.Hx, :, :] = nut.data[g.Nx+1:g.Nx+g.Hx, :, :] # west
        @views @. nut.data[:, 1:g.Hy, :] = nut.data[:, g.Ny+1:g.Ny+g.Hy, :] # south
        @views @. nut.data[:, :, 1:g.Hz] = nut.data[:, :, g.Hz+1:g.Hz+1] # bottom

        @views @. nut.data[g.Nx+g.Hx+1:g.Nx+2*g.Hx, :, :] = nut.data[1+g.Hx:2*g.Hx, :, :] # east
        @views @. nut.data[:, g.Ny+g.Hy+1:g.Ny+2*g.Hy, :] = nut.data[:, 1+g.Hy:2*g.Hy, :] # north
        @views @. nut.data[:, :, g.Nz+g.Hz+1:g.Nz+2*g.Hz] = nut.data[:, :, g.Nz+g.Hz:g.Nz+g.Hz] # top
    end
    return nothing
end
@inline function fill_halo_vel!(vels::NamedTuple, g::Grids)
    for vel in vels
        @views @. vel.data[1:g.Hx, :, :] = vel.data[g.Nx+1:g.Nx+g.Hx, :, :] # west
        @views @. vel.data[:, 1:g.Hy, :] = vel.data[:, g.Ny+1:g.Ny+g.Hy, :] # south
        @views @. vel.data[:, :, 1:g.Hz] = vel.data[:, :, g.Hz+1:g.Hz+1] # bottom

        @views @. vel.data[g.Nx+g.Hx+1:g.Nx+2*g.Hx, :, :] = vel.data[1+g.Hx:2*g.Hx, :, :] # east
        @views @. vel.data[:, g.Ny+g.Hy+1:g.Ny+2*g.Hy, :] = vel.data[:, 1+g.Hy:2*g.Hy, :] # north
        @views @. vel.data[:, :, g.Nz+g.Hz+1:g.Nz+2*g.Hz] = vel.data[:, :, g.Nz+g.Hz:g.Nz+g.Hz] * 0.4 # top
    end
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