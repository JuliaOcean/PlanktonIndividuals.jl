"""
    vel_copy!(vel::NamedTuple, u, v, w, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
Copy external velocities into `PlanktonModel`
"""
function vel_copy!(vel::NamedTuple, u, v, w, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    copy_interior_u!(vel.u.data, u, g, TX())
    copy_interior_v!(vel.v.data, v, g, TY())
    copy_interior_w!(vel.w.data, w, g, TZ())

    fill_halo_vel!(vel, g)
end

@inline function copy_interior!(c, t, g::AbstractGrid)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_u!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_v!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_w!(c, t, g::AbstractGrid, ::Periodic)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end

@inline function copy_interior_u!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx+1, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_v!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny+1, g.Hz+1:g.Hz+g.Nz), t)
end
@inline function copy_interior_w!(c, t, g::AbstractGrid, ::Bounded)
    copyto!(view(c, g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz+1), t)
end

function validate_temp(sim::PlanktonSimulation, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    temp_size = (g.Nx, g.Ny, g.Nz)
    validation = true

    if size(sim.input.temp)[1:3] ≠ temp_size
        validation = false
        throw(ArgumentError("Dimension mismatch: the size of temperature must be $(temp_size)."))
    end

    if floor(Int, sim.input.ΔT_temp / sim.ΔT) * sim.ΔT == sim.input.ΔT_temp
        nothing
    else
        throw(ArgumentError("ΔT_temp must be an integer times of $(sim.ΔT)."))
    end

    if size(sim.input.temp)[4] < floor(Int,sim.iterations*sim.ΔT/sim.input.ΔT_temp)
        throw(ArgumentError("Temperature provided cannot cover the time period of (iterations+1)*ΔT ($(sim.iterations*sim.ΔT) seconds) with ΔT_temp = $(sim.input.ΔT_temp)."))
    end

    return validation
end

function validate_PARF(sim::PlanktonSimulation, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    PARF_size = (g.Nx, g.Ny)
    validation = true

    if size(sim.input.PARF)[1:2] ≠ PARF_size
        validation = false
        throw(ArgumentError("Dimension mismatch: the size of PARF must be $(PARF_size)."))
    end

    if floor(Int, sim.input.ΔT_PAR / sim.ΔT) * sim.ΔT == sim.input.ΔT_PAR
        nothing
    else
        throw(ArgumentError("ΔT_PAR must be an integer times of $(sim.ΔT)."))
    end

    if size(sim.input.PARF)[3] < floor(Int,sim.iterations*sim.ΔT/sim.input.ΔT_PAR)
        throw(ArgumentError("Surface PAR provided cannot cover the time period of (iterations+1)*ΔT ($(sim.iterations*sim.ΔT) seconds) with ΔT_PAR = $(sim.input.ΔT_PAR)."))
    end

    return validation
end

function validate_velocity(sim::PlanktonSimulation, g::AbstractGrid{TX, TY, TZ}) where {TX, TY, TZ}
    if sim.input.vels ≠ (;)
        u_size = (g.Nx, g.Ny, g.Nz)
        v_size = (g.Nx, g.Ny, g.Nz)
        w_size = (g.Nx, g.Ny, g.Nz)
        validation = true

        if isa(TX(), Bounded)
            u_size = (g.Nx+1, g.Ny, g.Nz)
        end

        if isa(TY(), Bounded)
            v_size = (g.Nx, g.Ny+1, g.Nz)
        end

        if isa(TZ(), Bounded)
            w_size = (g.Nx, g.Ny, g.Nz+1)
        end

        if size(sim.input.vels.u)[1:3] ≠ u_size
            validation = false
            throw(ArgumentError("Dimension mismatch: the size of u must be $(u_size), for $(TX) topology."))
        end

        if size(sim.input.vels.v)[1:3] ≠ v_size
            validation = false
            throw(ArgumentError("Dimension mismatch: the size of v must be $(v_size), for $(TY) topology."))
        end

        if size(sim.input.vels.w)[1:3] ≠ w_size
            validation = false
            throw(ArgumentError("Dimension mismatch: the size of w must be $(w_size), for $(TZ) topology."))
        end

        if floor(Int, sim.input.ΔT_vel / sim.ΔT) * sim.ΔT == sim.input.ΔT_vel
            nothing
        else
            throw(ArgumentError("ΔT_vel must be an integer times of $(sim.ΔT)."))
        end

        if size(sim.input.vels.u)[4] < floor(Int,sim.iterations*sim.ΔT/sim.input.ΔT_vel)
            throw(ArgumentError("Velocities provided cannot cover the time period of (iterations+1)*ΔT ($(sim.iterations*sim.ΔT) seconds) with ΔT_vel = $(sim.input.ΔT_vel)."))
        end
        return validation
    else
        return true
    end
end

function set_vels_fields!(sim::PlanktonSimulation, uv, vv, wv) 
    sim.input.vels = (u = uv, v = vv, w = wv)
    validate_velocity(sim, sim.model.grid)
    return nothing
end
function set_PARF_fields!(sim::PlanktonSimulation, PARF) 
    sim.input.PARF = PARF
    validate_PARF(sim, sim.model.grid)
    return nothing
end
function set_temp_fields!(sim::PlanktonSimulation, temp) 
    sim.input.temp = temp
    validate_temp(sim, sim.model.grid)
    return nothing
end
