"""
update!(sim::PlanktonSimulaiton)
update the `PlanktonSimulaiton` for `sim.iterations` time steps.
"""
function update!(sim::PlanktonSimulation)
if sim.input.vels ≠ (;)
    for t in 1:sim.iterations
        t_vel = Int(sim.model.t÷sim.input.ΔT_vel)+1 # starting from 1
        vel_copy!(sim.model.timestepper.vel₁, sim.input.vels.u[:,:,:,t_vel],
                  sim.input.vels.v[:,:,:,t_vel], sim.input.vels.w[:,:,:,t_vel], sim.model.grid)

        # tc = sim.model.t % 86400 ÷ sim.ΔT + 1 # time of day in ΔT as Int index for temp and PAR
        t_par = Int(sim.model.t÷sim.input.ΔT_PAR)+1 # starting from 1
        copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,t_par])

        t_temp = Int(sim.model.t÷sim.input.ΔT_temp)+1 # starting from 1
        copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,t_temp], sim.model.grid)

        TimeStep!(sim.model, sim.ΔT, sim.diags)

        write_output!(sim.output_writer, sim.model, sim.diags, sim.ΔT)
    end
else
    for t in 1:sim.iterations
        t_par = Int(sim.model.t÷sim.input.ΔT_PAR)+1 # starting from 1
        copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,t_par])

        t_temp = Int(sim.model.t÷sim.input.ΔT_temp)+1 # starting from 1
        copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,t_temp], sim.model.grid)

        TimeStep!(sim.model, sim.ΔT, sim.diags)

        write_output!(sim.output_writer, sim.model, sim.diags, sim.ΔT)
    end
end
end