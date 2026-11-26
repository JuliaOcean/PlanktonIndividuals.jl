# #  Global Ocean Example
#
# Here we simulate phytoplankton cells as **passive** Lagrangian particles and nutrient fields as **passive** tracers in the global ocean.

# ## 1. Import packages

using PlanktonIndividuals, Plots, JLD2, PlanktonKernels

# ## 2. Generate Flow Fields
#
# First, we load velocity fields and a mask for the land shape information.
# Then we re-format velocity fields so they can be loaded in to PlanktonSimulation.

begin
    file = jldopen(PlanktonIndividuals.global_vels,"r")
    mask = Array(file["mask"] .!= 0)
    uvels = file["u"]
    vvels = file["v"]
    wvels = file["w"]
    close(file)
end

# Then we generate the grid information.

grid = LatLonGrid(size=(360, 160, 1), lat=(-80,80), lon=(-180,180), z=(0,-10), landmask = mask)

# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid, number of individuals, parameters, and nutrient initial conditions.
#
# In order to change individuals and tracer to **passive** mode, we need to update some parameters.

phyt_parameter = Dict("PCmax"   => [0.0], # maximum photosynthesis rate
                      "dvid_P"  => [0.0], # probability of cell division
                      "mort_P"  => [0.0], # probability of cell natural death
                      "thermal" => [0.0]  # thermal damage and repair
                     )

bgc_parameter = Dict("kDOC" => 0.0, # Remineralization rate for DOC
                     "Nit"  => 0.0, # Nitrification rate for NH4
                     "kDON" => 0.0, # Remineralization rate for DON
                     "kDOP" => 0.0, # Remineralization rate for DON
                     "kPOC" => 0.0, # Remineralization rate for POC
                     "kPON" => 0.0, # Remineralization rate for PON
                     "kPOP" => 0.0, # Remineralization rate for POP
                     "w_sink_org" => 0.0,
                     "w_sink_inorg" => 0.0,
                    )

# Then we update new parameter values in the model.

model = PlanktonModel(CPU(), grid;
                      mode = CarbonMode(),
                      phyto = phyto_setup(phyt_parameter, [360], 1),
                      options = model_options(max_individuals = 360*8),
                      bgc_params = bgc_parameter)

# We also need to setup a runtime simulation to run the model.
# The simulation includes time step, number of time steps, flow fields that
# will be used etc.

sim = PlanktonSimulation(model; ΔT = 3600.0,
                                iterations = 24,
                                vels=(u=uvels, v=vvels, w=wvels),
                                ΔT_vel=3600.0*24)

# ## 4. Run the Model
#
# Finaly, we run the model and plot the distribution of individuals as well as nutrient fields.
#
# We use Plots.jl to plot individuals and nutrient fields.

function plot_model(model::PlanktonModel, uu)
    ## Coordinate arrays for plotting
    xC, yC = collect(model.grid.xC)[3:end-2], collect(model.grid.yC)[3:end-2]

    ## heatmap of the flow field
    fl_plot = Plots.contourf(xC, yC, uu', xlabel="x (m)", ylabel="y (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 .- 180 # convert fractional indices to degree
    py = Array(model.individuals.phytos.sp1.data.y) .* 1 .- 80  # convert fractional indices to degree
    Plots.scatter!(fl_plot, px, py, ms=3, color = :red, legend=:none)

    ## DOC field
    trac1 = Plots.heatmap(xC, yC, Array(model.tracers.DOC.data)[3:end-2,3:end-2,3]', xlabel="x (m)", ylabel="y (m)", clims=(0.5, 1.1), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(1200, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "Passive Tracer"])

    return plt
end

# We run the model for 24 time steps (1 hour per time step) and plot the individuals and DOC field.

simulated_model = begin
    sim_run!(sim)
    sim.model
end

# We plot the current state of the model.

begin
    u_plot = uvels[:,:,1,1]
    u_plot[.!mask[:, :, 1]] .= NaN
    plot_model(simulated_model, u_plot)
end

# To animate the same 24-hour simulation, recreate the model and simulation,
# then use this loop in place of the model-run cell above.
#
# ```julia
# sim.iterations = 1
# anim = @animate for i in 1:24
#     sim_run!(sim)
#     plot_model(sim.model, uvels[:, :, 1, 1])
# end
# gif(anim; fps = 15)
# ```
# ![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_global.gif)
