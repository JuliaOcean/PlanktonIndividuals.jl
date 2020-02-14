using Random, Printf, Plots
using Oceananigans, Oceananigans.Utils
using PhytoAgentModel

### Oceananigans Setup ###
Nz = 32       # Number of grid points in x, y, z
Δz = 1.0      # Grid spacing in x, y, z (meters)
Qᵀ = 5e-5     # Temperature flux at surface
Qᵘ = -2e-5    # Velocity flux at surface
∂T∂z = 0.005    # Initial vertical temperature gradient
evaporation = 1e-7     # Mass-specific evaporation rate [m s⁻¹]
f = 1e-4     # Coriolis parameter
α = 2e-4     # Thermal expansion coefficient
β = 8e-4     # Haline contraction coefficient

u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ))

T_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Flux, Qᵀ),
                                bottom = BoundaryCondition(Gradient, ∂T∂z))

## Salinity flux: Qˢ = - E * S
@inline Qˢ(i, j, grid, time, iter, U, C, p) = @inbounds -p.evaporation * C.S[i, j, 1]

S_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qˢ))

model = Model(
         architecture = CPU(),
                 grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Δz*Nz, Δz*Nz, Δz*Nz)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = AnisotropicMinimumDissipation(),
  boundary_conditions = HorizontallyPeriodicSolutionBCs(u=u_bcs, T=T_bcs, S=S_bcs),
           parameters = (evaporation = evaporation,)
)

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density tradient with random noise superposed.
T₀(x, y, z) = 20 + ∂T∂z * z + ∂T∂z * model.grid.Lz * 1e-6 * Ξ(z)

## Velocity initial condition: random noise scaled by the friction velocity.
u₀(x, y, z) = sqrt(abs(Qᵘ)) * 1e-1 * Ξ(z)

set!(model, u=u₀, w=u₀, T=T₀, S=35)

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=5.0)

### warm up ###
for i in 1:20
    update_Δt!(wizard, model)

    ## Time step the model forward
    walltime = @elapsed time_step!(model, 20, wizard.Δt)

    ## Print a progress message
    @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
            model.clock.iteration, prettytime(model.clock.time), prettytime(wizard.Δt),
            prettytime(walltime))
end

### PlanktonAgents Setup ###
phy_grid = read_Ogrids(model.grid);

#                   Dim output, NutOutput, GridChoice, Gridoff, VelChoice, Veloff, SaveGrid, SaveVel, Test
RunOption=RunOptions(3, true,   true,      false,      Dict(),  false,     Dict(), false,    false,   false);
#                   Nindivi, Nsp, Nsuper,    Cquota(mmol/cell),  mean, var
PhytoOpt = PlankOpt(1000,    2,   Int(1e0),  [1.8e-11, 1.8e-10], 1.0,  0.25)

#                 Nindivi, Nsp, Nsuper,    Cquota(mmol/cell), mean, var
ZooOpt = PlankOpt(1000,    1,   Int(1e0),  [1.8e-9],          1.0,  0.25)

#                  nTime,  ΔT, PhytoOpt, Zoo,   ZooOpt
RunParam=RunParams(25*60,  60, PhytoOpt, true,  ZooOpt)

phy_model = PA_Model(phy_grid, RunParam;                  #DIC, DIN,  DOC,  DON, POC, PON, mmol/m3
                     nutrients = setup_nutrients(phy_grid,[2.0, 0.05, 20.0, 0.0, 0.0, 0.0]));

### run PlanktonAgents with velocities from Oceananigans ###
for i in 1:220
    vel_field =[]
    for j in 1:3
        u = model.velocities.u.data.parent
        v = model.velocities.v.data.parent
        w = model.velocities.w.data.parent
        vel = PhytoAgentModel.velocity(u, v, w)
        push!(vel_field,vel)
        time_step!(model, 6, 5)
    end
    PA_advectRK4!(phy_model, RunParam.ΔT, vel_field)
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
for i in 1:size(phy_model.individuals,1)
    convert_coordinates(B1[i],phy_grid)
    convert_coordinates(B2[i],phy_grid)
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
anim = @animate for i in 1:220
    p1 = scatter(B1[i].x,B1[i].y, zcolor=B1[i].size, m=(:diamond, 3, :algae, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), clims=(0,3), cbar=true, label = "sp1")
    p2 = scatter(B2[i].x,B2[i].y, zcolor=B2[i].size, m=(:circle, 3, :amp, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-1,33), clims=(0,3), cbar=true, label = "sp2")
    plt = plot(p1,p2,layout=grid(1,2),size=(800,300),dpi=100)
end
gif(anim,"tmp_xy2d.gif", fps = 15)

anim = @animate for i in 1:220
    p1 = scatter(B1[i].y, B1[i].z, zcolor=B1[i].size, m=(:diamond, 3, :algae, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-50,5), cbar = true, clims=(0,3), label = "sp1")
    p2 = scatter(B2[i].y, B2[i].z, zcolor=B2[i].size, m=(:circle, 3, :amp, 0.8, Plots.stroke(0)), xlims=(-1,33), ylims=(-50,5), clims=(0,3), cbar = true, label = "sp2")
    plt = plot(p1,p2,layout=grid(1,2),size=(800,300),dpi=100)
end
gif(anim,"tmp_yz2d.gif", fps = 15)
