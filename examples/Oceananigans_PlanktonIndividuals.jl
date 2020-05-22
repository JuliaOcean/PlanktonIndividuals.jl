using Random, Printf
using Oceananigans, Oceananigans.Utils
using PlanktonIndividuals

Nx = 25       # Number of grid points in x and y
Nz = 50       # Number of grid points in z
Δz = 4.0      # Grid spacing in x, y, z (meters)
Qᵀ = 1e-5     # Temperature flux at surface
∂T∂z = 0.005  # Initial vertical temperature gradient
f = 1e-4      # Coriolis parameter
α = 2e-4      # Thermal expansion coefficient
β = 8e-4      # Haline contraction coefficient

grid = RegularCartesianGrid(size=(Nx, Nx, Nz), extent=(Δz*Nx, Δz*Nx, Δz*Nz))
T_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵀ), bottom = BoundaryCondition(Gradient, ∂T∂z))

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

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=60.0);

resultpath = PrepRunDir()
phy_grid = read_Ogrids(model.grid, resultpath);
RunParam.nTime = 1440*10
RunParam.ΔT = 60
#           DIC  NH4   NO3   PO4   DOC  DON  DOP  POC  PON  POP  ZOO
nut_init = [2.0, 0.50, 1.20, 0.20, 1.0, 0.5, 0.1, 1.0, 0.5, 0.1, 2.0];

# add diagnostics
RunParam.params["diag_inds"] = [1,1,1,1,1,1,1,1,1,1, 1, 1, 1, 1]

phy_model = PI_Model(phy_grid, RunParam;
                     nutrients = setup_nutrients(phy_grid, nut_init));

simulation = Simulation(model, Δt=wizard, stop_iteration=40*20, progress_frequency=20)
run!(simulation)

Nsimulation = Simulation(model, Δt=10.0, stop_iteration=802, progress_frequency=2)

for i in 1:RunParam.nTime
    vel_field =[]
    for j in 1:3
        run!(Nsimulation)
        Nsimulation.stop_iteration += 2
        u = model.velocities.u.data.parent
        v = model.velocities.v.data.parent
        w = model.velocities.w.data.parent
        vel = PlanktonIndividuals.velocity(u, v, w)
        push!(vel_field,vel)
    end
    vel_itps = (generate_vel_itp(phy_model.grid, vel_field[1]),
                generate_vel_itp(phy_model.grid, vel_field[2]),
                generate_vel_itp(phy_model.grid, vel_field[3]))
    PI_advectRK4!(phy_model, RunParam.ΔT, vel_itps)
    PI_TimeStep!(phy_model, RunParam.ΔT, vel_field[end], resultpath)

    phyts_sp = sort_species(phy_model.individuals.phytos, phy_model.params["P_Nsp"])
    write_species_dynamics(phy_model.t, phyts_sp, resultpath)
    # write_output(phyts_sp, resultpath, phy_model.t)
    # write_nut_nc_each_step(phy_model.nutrients, phy_model.t, resultpath)
end

## Plots
using PyPlots, DelimitedFiles

## temporal variations of population and nutrient pools of individuals
Ndata = readdlm("results/cons_N.txt")
Pdata = readdlm("results/cons_P.txt")
sp1_data = readdlm("results/dynamic_species001.txt")

fig, axs = PyPlot.subplots(4,2,sharex=true,figsize=(12,10))
axs[1].plot(sp1_data[:,2],label ="Population")
axs[2].plot(sp1_data[:,5],label ="Average Size")
axs[3].plot(sp1_data[:,3],label ="Generation")
axs[4].plot(sp1_data[:,6],label ="C Biomass")
axs[4].plot(sp1_data[:,7],label ="C Reserve")
axs[5].plot(sp1_data[:,8],label ="N Reserve")
axs[6].plot(sp1_data[:,9],label ="P Reserve")
axs[7].plot(Ndata[:,4],label ="NH4")
axs[7].plot(Ndata[:,5],label ="NO3")
axs[8].plot(Ndata[:,4],label ="PO4")

## vertical and temporal variations of division rate
pops[:,:,1]= sum(sum(phy_model.diags.pop[:,:,:,:,1,1],dims=1),dims=2)[1,1,:,:]
pops[:,:,2]= sum(sum(phy_model.diags.tr[:,:,:,:,2] ./ 60,dims=1),dims=2)[1,1,:,:]
dvi_ver = mean(pops[:,:,1],dims=2)[:,1]
dvi_t = mean(pops[:,:,1],dims=1)[1,:]
pop_ver = mean(pops[:,:,2],dims=2)[:,1]
pop_t = mean(pops[:,:,2],dims=1)[1,:]

## vertical profile of division rate
fig,ax = plt.subplots(figsize=(4,7))
ax.plot(dvi_ver ./ pop_ver .* 24,collect(-199:4:-1),label = "Species 1")
ax.legend(loc=2, fontsize=8)
ax.set_xlabel("Division (per day)")
ax.set_ylabel("Depth (m)");

## temporal variations of division rate
fig,ax = plt.subplots(figsize=(7,4))
ax.plot(collect(1:1:240),dvi_t ./ pop_t .* 24,label = "Species 1")
ax.legend(loc=2, fontsize=10)
ax.set_ylabel("Division (per day)")
ax.set_xlabel("Time (Day)");
