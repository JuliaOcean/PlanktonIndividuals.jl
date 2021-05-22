##### fill halo points based on topology
@inline  fill_halo_west!(c, H::Int64, N::Int64, ::Periodic)  = @views @. c[1:H, :, :] = c[N+1:N+H, :, :]
@inline fill_halo_south!(c, H::Int64, N::Int64, ::Periodic)  = @views @. c[:, 1:H, :] = c[:, N+1:N+H, :]
@inline   fill_halo_top!(c, H::Int64, N::Int64, ::Periodic)  = @views @. c[:, :, 1:H] = c[:, :, N+1:N+H]

@inline   fill_halo_east!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[N+H+1:N+2H, :, :] = c[1+H:2H, :, :]
@inline  fill_halo_north!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[:, N+H+1:N+2H, :] = c[:, 1+H:2H, :]
@inline fill_halo_bottom!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[:, :, N+H+1:N+2H] = c[:, :, 1+H:2H]

@inline  fill_halo_west!(c, H::Int64, N::Int64, ::Bounded)   = @views @. c[1:H, :, :] = c[H+1:H+1, :, :]
@inline fill_halo_south!(c, H::Int64, N::Int64, ::Bounded)   = @views @. c[:, 1:H, :] = c[:, H+1:H+1, :]
@inline   fill_halo_top!(c, H::Int64, N::Int64, ::Bounded)   = @views @. c[:, :, 1:H] = c[:, :, H+1:H+1]

@inline   fill_halo_east!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[N+H+1:N+2H, :, :] = c[N+H:N+H, :, :]
@inline  fill_halo_north!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[:, N+H+1:N+2H, :] = c[:, N+H:N+H, :]
@inline fill_halo_bottom!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[:, :, N+H+1:N+2H] = c[:, :, N+H:N+H]

@inline   fill_halo_east_vel!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[N+H+2:N+2H, :, :] = c[N+H+1:N+H+1, :, :]
@inline  fill_halo_north_vel!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[:, N+H+2:N+2H, :] = c[:, N+H+1:N+H+1, :]
@inline fill_halo_bottom_vel!(c, H::Int64, N::Int64, ::Bounded)  = @views @. c[:, :, N+H+2:N+2H] = c[:, :, N+H+1:N+H+1]

@inline   fill_halo_east_Gc!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[N+H+1:N+2H, :, :] = c[1+H:2H, :, :]
@inline  fill_halo_north_Gc!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[:, N+H+1:N+2H, :] = c[:, 1+H:2H, :]
@inline fill_halo_bottom_Gc!(c, H::Int64, N::Int64, ::Periodic) = @views @. c[:, :, N+H+1:N+2H] = c[:, :, 1+H:2H]

@inline   fill_halo_east_Gc!(c, H::Int64, N::Int64, ::Bounded) = @views @. c[N+H+1:N+2H, :, :] = 0.0
@inline  fill_halo_north_Gc!(c, H::Int64, N::Int64, ::Bounded) = @views @. c[:, N+H+1:N+2H, :] = 0.0
@inline fill_halo_bottom_Gc!(c, H::Int64, N::Int64, ::Bounded) = @views @. c[:, :, N+H+1:N+2H] = 0.0

@inline function fill_halo_nut!(nuts::NamedTuple, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    for nut in nuts
          fill_halo_west!(nut.data, g.Hx, g.Nx, TX())
          fill_halo_east!(nut.data, g.Hx, g.Nx, TX())
         fill_halo_south!(nut.data, g.Hy, g.Ny, TY())
         fill_halo_north!(nut.data, g.Hy, g.Ny, TY())
           fill_halo_top!(nut.data, g.Hz, g.Nz, TZ())
        fill_halo_bottom!(nut.data, g.Hz, g.Nz, TZ())
    end
    return nothing
end

@inline function fill_halo_Gcs!(nuts::NamedTuple, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    for nut in nuts
          fill_halo_east_Gc!(nut.data, g.Hx, g.Nx, TX())
         fill_halo_north_Gc!(nut.data, g.Hy, g.Ny, TY())
        fill_halo_bottom_Gc!(nut.data, g.Hz, g.Nz, TZ())
    end
    return nothing
end

@inline function fill_halo_u!(u, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    if isa(TX(), Bounded)
        fill_halo_east_vel!(u, g.Hx, g.Nx, TX())
    else
        fill_halo_east!(u, g.Hx, g.Nx, TX())
    end
      fill_halo_west!(u, g.Hx, g.Nx, TX())
     fill_halo_south!(u, g.Hy, g.Ny, TY())
     fill_halo_north!(u, g.Hy, g.Ny, TY())
       fill_halo_top!(u, g.Hz, g.Nz, TZ())
    fill_halo_bottom!(u, g.Hz, g.Nz, TZ())
end

@inline function fill_halo_v!(v, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    if isa(TY(), Bounded)
        fill_halo_north_vel!(v, g.Hy, g.Ny, TY())
    else
        fill_halo_north!(v, g.Hy, g.Ny, TY())
    end
      fill_halo_west!(v, g.Hx, g.Nx, TX())
      fill_halo_east!(v, g.Hx, g.Nx, TX())
     fill_halo_south!(v, g.Hy, g.Ny, TY())
       fill_halo_top!(v, g.Hz, g.Nz, TZ())
    fill_halo_bottom!(v, g.Hz, g.Nz, TZ())
end

@inline function fill_halo_w!(w, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    if isa(TZ(), Bounded)
        fill_halo_bottom_vel!(w, g.Hz, g.Nz, TZ())
    else
        fill_halo_bottom!(w, g.Hz, g.Nz, TZ())
    end
     fill_halo_west!(w, g.Hx, g.Nx, TX())
     fill_halo_east!(w, g.Hx, g.Nx, TX())
    fill_halo_south!(w, g.Hy, g.Ny, TY())
    fill_halo_north!(w, g.Hy, g.Ny, TY())
      fill_halo_top!(w, g.Hz, g.Nz, TZ())
end

@inline function fill_halo_vel!(vels::NamedTuple, g::AbstractGrid)
    fill_halo_u!(vels.u.data, g)
    fill_halo_v!(vels.v.data, g)
    fill_halo_w!(vels.w.data, g)
    return nothing
end
