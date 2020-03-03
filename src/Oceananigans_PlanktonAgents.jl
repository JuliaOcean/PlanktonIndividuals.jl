using Random, Printf, Plots
using Oceananigans, Oceananigans.Utils
using PhytoAgentModel

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
                 grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Δz*Nz, Δz*Nz, Δz*Nz)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = AnisotropicMinimumDissipation(),
  boundary_conditions = (T=T_bcs,),
           parameters = (evaporation = evaporation,)
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

#                   Dim output, NutOutput, GridChoice, Gridoff, VelChoice, Veloff, SaveGrid, SaveVel, Test
RunOption=RunOptions(3, true,   true,      false,      Dict(),  false,     Dict(), false,    false,   false);
#                   Nindivi, Nsp, Nsuper,    Cquota(mmol/cell),  mean, var
PhytoOpt = PlankOpt(1000,    2,   Int(1e0),  [1.8e-11, 1.8e-10], 1.0,  0.25)

#                 Nindivi, Nsp, Nsuper,    Cquota(mmol/cell), mean, var
#ZooOpt = PlankOpt(1000,    1,   Int(1e0),  [1.8e-9],          1.0,  0.25)

#                  nTime,  ΔT, PhytoOpt, Zoo,   ZooOpt
RunParam=RunParams(25*60,  60, PhytoOpt, false, nothing)

#           DIC, DIN,  DOC,  DON, POC, PON, mmol/m3
nut_init = [2.0, 0.05, 20.0, 0.0, 0.0, 0.0];
phy_model = PA_Model(phy_grid, RunParam; nutrients = setup_nutrients(phy_grid, nut_init));

### run PlanktonAgents with velocities from Oceananigans ###
Nsimulation = Simulation(model, Δt=5.0, stop_iteration=506, progress_frequency=6)
for i in 1:500
    vel_field =[]
    for j in 1:3
        u = model.velocities.u.data.parent
        v = model.velocities.v.data.parent
        w = model.velocities.w.data.parent
        vel = PhytoAgentModel.velocity(u, v, w)
        push!(vel_field,vel)
        run!(Nsimulation)
        Nsimulation.stop_iteration += 6
    end
    vel_itps = (generate_vel_itp(model.grid, vel_field[1]),
                generate_vel_itp(model.grid, vel_field[2]),
                generate_vel_itp(model.grid, vel_field[3]))
    PA_advectRK4!(phy_model, RunParam.ΔT, vel_itps)
    PA_TimeStep!(phy_model, RunParam.ΔT, vel_field[end])
end

### Post-processing ###
B1 = []; B2 = [];
for i in 1:size(phy_model.individuals,1)
    sort_species(phy_model.individuals[i,1], B1, 1)
    sort_species(phy_model.individuals[i,1], B2, 2)
end
HD1 = []; HD2 = [];
for i in 1:size(phy_model.individuals,1)
    HD_1 = count_horizontal_num(B1[i],phy_grid);
    push!(HD1,HD_1)
    HD_2 = count_horizontal_num(B2[i],phy_grid);
    push!(HD2,HD_2)
end
output1 = compute_mean_species(B1, size(phy_model.individuals,1));
output2 = compute_mean_species(B2, size(phy_model.individuals,1));

VD1 = []; VD2 = [];
for i in 1:size(phy_model.individuals,1)
    VD_1 = count_vertical_num(B1[i]);
    push!(VD1,VD_1)
    VD_2 = count_vertical_num(B2[i]);
    push!(VD2,VD_2)
end

### Plots & Animations ###
anim = @animate for i in 1:500
    p1 = scatter(B1[i].x,B1[i].y,B1[i].z, zcolor=B1[i].size, m=(:diamond, 3, :algae, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), zlims=(-50,5), clims=(0,3), cbar=true, label = "species1")
    scatter!(p1,Zoos[i].x,Zoos[i].y,Zoos[i].z, zcolor=Zoos[i].size, m=(:star8, 5, :ice, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), zlims=(-50,5), clims=(0,3), cbar=false, label = "zoo")
    scatter!(p1,B2[i].x,B2[i].y,B2[i].z, zcolor=B2[i].size, m=(:circle, 3, :amp, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), zlims=(-50,5), clims=(0,3), cbar=false, label = "species2")
    plt = plot(p1,size=(800,500),dpi=100)
end
gif(anim,"tmp_test.gif", fps = 15)

anim = @animate for i in 1:500
    p1 = scatter(B1[i].x,B1[i].y, zcolor=B1[i].size, m=(:diamond, 3, :algae, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), clims=(0,3), cbar=true, label = "sp1")
    p2 = scatter(B2[i].x,B2[i].y, zcolor=B2[i].size, m=(:circle, 3, :amp, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), clims=(0,3), cbar=true, label = "sp2")
    plt = plot(p1,p2,layout=grid(1,2),size=(800,300),dpi=100)
end
gif(anim,"tmp_xy2d.gif", fps = 15)

anim = @animate for i in 1:500
    p1 = scatter(B1[i].y, B1[i].z, zcolor=B1[i].size, m=(:diamond, 3, :algae, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-50,5), cbar = true, clims=(0,3), label = "sp1")
    p2 = scatter(B2[i].y, B2[i].z, zcolor=B2[i].size, m=(:circle, 3, :amp, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-50,5), clims=(0,3), cbar = true, label = "sp2")
    plt = plot(p1,p2,layout=grid(1,2),size=(800,300),dpi=100)
end
gif(anim,"tmp_yz2d.gif", fps = 15)
