mutable struct PI_Input
    temp::AbstractArray{Float64,4}      # temperature
    PARF::AbstractArray{Float64,3}      # PARF
    vels::NamedTuple                    # velocity fields for nutrients and individuals
end

mutable struct PI_simulation
    model::PI_Model                  # Model object
    input::PI_Input                  # model input, temp, PAR, and velocities
    ΔT::Int64                        # model time step
    nΔT::Int64                       # number of time steps to run in the simulation
    diag_freq::Int64                 # frequency of diagnostics
    res_dir::Union{String,Nothing}   # directory to store results
    save_diags::Bool
    save_individuals::Bool
    vel_reuse::Bool
end

"""
    PI_simulation(model; ΔT, nΔT, diag_freq,
                  PARF_path = dirname(pathof(PlanktonIndividuals))*"/../samples/PAR.bin",
                  temp_path = dirname(pathof(PlanktonIndividuals))*"/../samples/temp.bin",
                  vels = (;),
                  res_dir = nothing
                  save_diags = false,
                  save_individuals = false,
                  vel_reuse = false
                  )
Generate the `PI_simulation` struct for a `model` with time step `ΔT`. 

Keyword Arguments
=================
- `ΔT` (required): `model` time step in second.
- `nΔT` (required): The number of time steps to run in `simulation`.
- `diag_freq` (required): The frequency of diagnostics (in number of time steps).
- `PARF_path` and `temp_path` (optional): External forcings of PAR and temperature.
- `res_dir` (optional): Create a directory to store results, `nothing` by default.
- `save_diags` and `save_individuals` (optional): whether to save diagnostics or individuals.
- `vels` (optional): The velocity fields for nutrient fields and individuals. `nothing` means no velocities
                     will be applied in the simulation. Otherwise, `vels` mush be a `NamedTuple` containing
                     all `u`, `v`, and `w`. Each of `u`, `v`, and `w` must be an 4D-`Array` of 
                     `(Nx+2Hx, Ny+2Hy, Nz+2Hz, nΔT+1)` elements, including halo points and initial conditions.
"""
function PI_simulation(model::PI_Model; ΔT::Int64, nΔT::Int64, diag_freq::Int64,
                       PARF_path = PlanktonIndividuals.default_PAR,
                       temp_path = PlanktonIndividuals.default_temperature,
                       vels = (;),
                       res_dir = nothing,
                       save_diags = false,
                       save_individuals = false,
                       vel_reuse = false
                       )

    if vels ≠ (;)
        grid_size = (model.grid.Nx, model.grid.Ny, model.grid.Nz)
        grid_size_w = (model.grid.Nx, model.grid.Ny, model.grid.Nz+1)
        if size(vels.u)[1:3] == size(vels.v)[1:3] == grid_size
            if size(vels.w)[1:3] == grid_size_w
                vel_copy!(model.timestepper.vel₀, vels.u[:,:,:,1],
                        vels.v[:,:,:,1], vels.w[:,:,:,1], model.grid)
            else
                throw(ArgumentError("Dimension mismatch: the size of w must be $(grid_size_w), Nz+1 layers for bounded direction."))
            end
        else
            throw(ArgumentError("Dimension mismatch: the size of u and v must be $(grid_size)."))
        end

        if size(vels.u)[4] < nΔT
            throw(ArgumentError("Velocities provided not enough for $(nΔT) time steps."))
        end

    end

    temp = read_temp_input(ΔT, model.grid, path = temp_path)
    PARF = read_IR_input(ΔT, model.grid, path = PARF_path)
    input = PI_Input(temp, PARF, vels)

    sim = PI_simulation(model, input, ΔT, nΔT, diag_freq, res_dir, save_diags, save_individuals, vel_reuse)

    return sim
end

import Base: show

function show(io::IO, sim::PI_simulation)
    print(io, "ΔT: $(sim.ΔT)s\n",
              "model time: $(sim.model.t)s\n",
              "number of time steps: $(sim.nΔT)\n",
              "save averaged diagnostics every $(sim.diag_freq*sim.ΔT)s, $(sim.diag_freq) time steps\n",
              "results saved at $(sim.res_dir)\n",
              "save diags: $(sim.save_diags)\n",
              "save individuals: $(sim.save_individuals)\n",
              "reuse velocity fields: $(sim.vel_reuse)\n")
end


"""
    update!(sim::PI_simulaiton)
update the `PI_simulaiton` for `sim.nΔT` time steps.
"""
function update!(sim::PI_simulation)
    if sim.input.vels ≠ (;)
        if sim.vel_reuse
            ti = 1
            te = sim.nΔT
        else
            ti = sim.model.t÷sim.ΔT + 1
            te = sim.model.t÷sim.ΔT + sim.nΔT
        end
        for t in ti:te
            tc = sim.model.t % 86400 ÷ sim.ΔT + 1 # time of day in ΔT as Int index for temp and PAR
            vel_copy!(sim.model.timestepper.vel₁, sim.input.vels.u[:,:,:,t],
                      sim.input.vels.v[:,:,:,t], sim.input.vels.w[:,:,:,t], sim.model.grid)
            copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,tc])
            copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,tc], sim.model.grid)

            if sim.res_dir ≠ nothing
                PI_TimeStep!(sim.model, sim.ΔT, sim.res_dir)
                if sim.save_diags
                    if sim.model.t % (sim.diag_freq*sim.ΔT) == 0.0
                        write_diags_to_jld2(sim.model.diags, sim.res_dir, sim.model.t, sim.diag_freq)
                    end
                end
                if sim.save_individuals
                    write_individuals_to_jld2(sim.model.individuals.phytos, sim.res_dir, sim.model.t,
                                              atts=(:x, :y, :z, :Sz))
                end
            else
                PI_TimeStep!(sim.model, sim.ΔT)
            end
        end
    else
        for t in 1:sim.nΔT
            tc = sim.model.t % 86400 ÷ sim.ΔT + 1 # time of day in ΔT as Int index for temp and PAR
            copyto!(sim.model.timestepper.PARF, sim.input.PARF[:,:,tc])
            copy_interior!(sim.model.timestepper.temp, sim.input.temp[:,:,:,tc], sim.model.grid)
            if sim.res_dir ≠ nothing
                PI_TimeStep!(sim.model, sim.ΔT, sim.res_dir)
                if sim.save_diags
                    if sim.model.t % (sim.diag_freq*sim.ΔT) == 0.0
                        write_diags_to_jld2(sim.model.diags, sim.res_dir, sim.model.t, sim.diag_freq)
                    end
                end
                if sim.save_individuals
                    write_individuals_to_jld2(sim.model.individuals.phytos, sim.res_dir, sim.model.t,
                                              atts=(:x, :y, :z, :Sz))
                end
            else
                PI_TimeStep!(sim.model, sim.ΔT)
            end
        end
    end
end
