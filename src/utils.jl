##### copy external velocities into the model
function vel_copy!(vel::NamedTuple, u, v, w, g::AbstractGrid)
    copy_interior!(vel.u.data, u, g)
    copy_interior!(vel.v.data, v, g)
    copy_interior_bounded!(vel.w.data, w, g)

    fill_halo_vel!(vel, g)
end

function zero_fields!(a)
    for i in 1:length(a)
        @inbounds a[i].data .= 0.0
    end
end

@inline function copy_interior!(c, t, g::AbstractGrid)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_bounded!(c, t, g::AbstractGrid)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz+1), t)
end


