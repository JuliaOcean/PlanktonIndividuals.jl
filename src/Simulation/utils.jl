function validate_temp(sim::PlanktonSimulation, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
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

function validate_PARF(sim::PlanktonSimulation, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
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

function validate_velocity(sim::PlanktonSimulation, g::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
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
    FT = sim.model.FT 
    sim.input.vels = (u = FT.(uv), v = FT.(vv), w = FT.(wv))
    validate_velocity(sim, sim.model.grid)
    return nothing
end
function set_PARF_fields!(sim::PlanktonSimulation, PARF) 
    FT = sim.model.FT 
    sim.input.PARF = FT.(PARF)
    validate_PARF(sim, sim.model.grid)
    return nothing
end
function set_temp_fields!(sim::PlanktonSimulation, temp) 
    FT = sim.model.FT 
    sim.input.temp = FT.(temp)
    validate_temp(sim, sim.model.grid)
    return nothing
end
