"""
    update!(sim::PlanktonSimulation; time_offset = (vels = true, PAFR = true, temp = true))

Update the `PlanktonSimulation` for `sim.iterations` time steps.
`time_offset` is used when velocities (or PARF or temperature) starts from timestep 1, but model.t is not.
It is usually used when velocity fields are too large and need to be broken down into several parts.
Only one part of the whole velocity fields can be constructed into a `PlanktonSimulation`, so in this PlanktonSimulation
`model.iteration` might no be 1, but the velocity fields need to start from 1 (same for PARF or temperature fields).

"""
function update!(sim::PlanktonSimulation; time_offset = (vels = false, PARF = false, temp = false))
    if sim.input.vels ≠ (;)
        for i in 1:sim.iterations
            model_t_vels = time_offset.vels ? (i-1)*sim.ΔT : sim.model.t
            model_t_PARF = time_offset.PARF ? (i-1)*sim.ΔT : sim.model.t
            model_t_temp = time_offset.temp ? (i-1)*sim.ΔT : sim.model.t

            t_vel = floor(Int, model_t_vels/sim.input.ΔT_vel)+1 # starting from 1
            vel_copy!(sim.model.timestepper.vel₁, sim.input.vels.u[:,:,:,t_vel],
                    sim.input.vels.v[:,:,:,t_vel], sim.input.vels.w[:,:,:,t_vel], sim.model.grid)

            t_par = floor(Int,model_t_PARF/sim.input.ΔT_PAR)+1 # starting from 1
            copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,t_par])

            t_temp = floor(Int,model_t_temp/sim.input.ΔT_temp)+1 # starting from 1
            copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,t_temp], sim.model.grid)

            TimeStep!(sim.model, sim.ΔT, sim.diags)

            write_output!(sim.output_writer, sim.model, sim.diags, sim.ΔT)
        end
    else
        for i in 1:sim.iterations
            model_t_PARF = time_offset.PARF ? sim.model.t : (i-1)*sim.ΔT
            model_t_temp = time_offset.temp ? sim.model.t : (i-1)*sim.ΔT

            t_par = floor(Int,model_t_PARF/sim.input.ΔT_PAR)+1 # starting from 1
            copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,t_par])

            t_temp = floor(Int,model_t_temp/sim.input.ΔT_temp)+1 # starting from 1
            copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,t_temp], sim.model.grid)

            TimeStep!(sim.model, sim.ΔT, sim.diags)

            write_output!(sim.output_writer, sim.model, sim.diags, sim.ΔT)
        end
    end
end