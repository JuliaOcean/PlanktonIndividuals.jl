using Random, Printf
using Oceananigans, Oceananigans.Utils
using PlanktonIndividuals

### Oceananigans Setup with heat induced convection ###
Nz = 32       # Number of grid points in x, y, z
Δz = 1.0      # Grid spacing in x, y, z (meters)
Qᵀ = 5e-5     # Temperature flux at surface
∂T∂z = 0.005    # Initial vertical temperature gradient
evaporation = 1e-7     # Mass-specific evaporation rate [m s⁻¹]
f = 1e-4     # Coriolis parameter
α = 2e-4     # Thermal expansion coefficient
β = 8e-4     # Haline contraction coefficient

grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Δz*Nz, Δz*Nz, Δz*Nz))
T_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵀ),
                                 bottom = BoundaryCondition(Gradient, ∂T∂z))

model = IncompressibleModel(
         architecture = CPU(),
                 grid = grid,
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = AnisotropicMinimumDissipation(),
  boundary_conditions = (T=T_bcs,)
)

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density tradient with random noise superposed.
T₀(x, y, z) = 20 + ∂T∂z * z + ∂T∂z * model.grid.Lz * 1e-6 * Ξ(z)

set!(model, T=T₀)

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=5.0)

### warm up ###
simulation = Simulation(model, Δt=wizard, stop_iteration=25*20, progress_frequency=20)
run!(simulation)

### PlanktonAgents Setup ###
phy_grid = read_Ogrids(model.grid);

#                   Nindivi, Nsp, Nsuper,    Cquota(mmol/cell),  mean, var
PhytoOpt = PlankOpt(3,    2,   Int(1e0),  [1.8e-11, 1.8e-10], 1.0,  0.25)

#                 Nindivi, Nsp, Nsuper,    Cquota(mmol/cell), mean, var
#ZooOpt = PlankOpt(1000,    1,   Int(1e0),  [1.8e-9],          1.0,  0.25)

#                  nTime,  ΔT, PhytoOpt, Zoo,   ZooOpt
RunParam=RunParams(25*60,  60, PhytoOpt, false, nothing)

#           DIC, NH4, NO3, PO4, DOC,  DON, PON, POC, PON, POP mmol/m3
nut_init = [2.0, 0.05,0.05,0.01,20.0, 0.0, 0.0, 0.0, 0.0, 0.0];
phy_model = PI_Model(phy_grid, RunParam;
                     nutrients = setup_nutrients(phy_grid, nut_init));
resultpath = PrepRunDir()

### run PlanktonAgents with velocities from Oceananigans ###
Nsimulation = Simulation(model, Δt=5.0, stop_iteration=506, progress_frequency=6)
for i in 1:500
    vel_field =[]
    for j in 1:3
        u = model.velocities.u.data.parent
        v = model.velocities.v.data.parent
        w = model.velocities.w.data.parent
        vel = PlanktonIndividuals.velocity(u, v, w)
        push!(vel_field,vel)
        run!(Nsimulation)
        Nsimulation.stop_iteration += 6
    end
    vel_itps = (generate_vel_itp(phy_model.grid, vel_field[1]),
                generate_vel_itp(phy_model.grid, vel_field[2]),
                generate_vel_itp(phy_model.grid, vel_field[3]))
    PI_advectRK4!(phy_model, RunParam.ΔT, vel_itps)
    PI_TimeStep!(phy_model, RunParam.ΔT, vel_field[end], resultpath)
    write_output(phy_model.individuals, resultpath, phy_model.t*RunParam.ΔT)
end
