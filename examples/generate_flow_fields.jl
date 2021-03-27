using Random, Statistics
using Oceananigans

Nx = 32       # Number of grid points in x 
Ny = 32       # Number of grid points in y
Nz = 32       # Number of grid points in z
Δz = 4.0      # Grid spacing in x, y, z (meters)

grid = RegularRectilinearGrid(size = (Nx, Ny, Nz), extent=(Δz*Nx, Δz*Ny, Δz*Nz),topology = (Periodic, Periodic, Bounded))

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=2e-4, β=8e-4))

Qʰ = 200  # W m⁻², surface _heat_ flux
ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991 # J K⁻¹ s⁻¹, typical heat capacity for seawater

Qᵀ = Qʰ / (ρₒ * cᴾ) # K m⁻¹ s⁻¹, surface _temperature_ flux

dTdz = 0.01 # K m⁻¹

T_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵀ),  bottom = BoundaryCondition(Gradient, dTdz))

u₁₀ = 1.00    # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3  # dimensionless drag coefficient
ρₐ = 1.225   # kg m⁻³, average density of air at sea-level

Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ))

@inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate * S

evaporation_rate = 1e-3 / 3600 

evaporation_bc = BoundaryCondition(Flux, Qˢ, field_dependencies=:S, parameters=evaporation_rate)

S_bcs = TracerBoundaryConditions(grid, top=evaporation_bc)

model = IncompressibleModel(architecture = CPU(),
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            coriolis = FPlane(f=1e-4),
                            buoyancy = buoyancy,
                            closure = AnisotropicMinimumDissipation(),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs))

Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)

uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=35)

wizard = TimeStepWizard(cfl=1.0, Δt=1.0, max_change=1.1, max_Δt=20.0)

simulation = Simulation(model, iteration_interval = 10, stop_time=60*60*2, Δt = wizard)
run!(simulation)

Nsimulation = Simulation(model, iteration_interval = 10, stop_time=60*60*3, Δt = 20.0)

Nsimulation.output_writers[:fields] =
        JLD2OutputWriter(model, model.velocities,
                         schedule = TimeInterval(60),
                         dir = ".",
                         prefix = "velocities",
                         force = true)

run!(Nsimulation)