# # Horizontal 2-Dimensional Example
# Here we simulate phytoplankton cells as Lagrangian particles in a horizontal, two-dimensional, flow field.
# The domain is periodic in both directions.

# ## 1. Import packages

using PlanktonKernels.Units
using PlanktonIndividuals, Plots, PlanktonKernels

include(joinpath(pkgdir(PlanktonIndividuals), "examples", "helper_functions.jl"))

# ## 2. Generate Flow Fields
# First we generate grid information and the computational architecture (CPU).

arch = CPU()

grid = RectilinearGrid(size=(128, 128, 1), x = (0.0,128meters), y = (0.0,128meters), z = (0.0,-1meter))

# Then we use a stream function to generate the flow field which is a double-gyre configuration [as explained here](https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html#Sec7.1)

(uvels, vvels, wvels, ϕcenters) = streamfunction_xy();

# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid, and plankton community.

model = PlanktonModel(arch, grid; phyto = phyto_setup(nothing, [2^7], 1),
                                  options = model_options(max_individuals = 2^7*8))

# We modify the initial tracer field to highlight advection.
#
# !!! note
#     DOC is used here as a generic tracer; not to represent dissolved organic carbon per se.

initialized_model = begin
    for i in 1:132
        for j in 1:132
            model.tracers.DOC.data[i,j,:] .= isless(132-j,i) .+ 1
        end
    end
    model
end

# Finally we setup the duration of the model simulation and the kind of output we want.

sim = PlanktonSimulation(initialized_model; ΔT = 60.0, iterations = 1,
                                vels=(u=uvels, v=vvels, w=wvels),
                                ΔT_vel=60.0*240)

# ## 4. Model Run
#
# We run the model for 240 time steps (4 hours) and generate animation.

# To plot the distribution of individuals as well as nutrient fields we use Plots.jl and create a function that can easily be re-used e.g. to create an animation.

function plot_model(model::PlanktonModel)
    ## coordinate arrays for plotting
    xC, yC = collect(model.grid.xC)[3:130], collect(model.grid.yC)[3:130]

    ## heatmap of the flow field
    fl_plot = Plots.contourf(xC, yC, ϕcenters', xlabel="x (m)", ylabel="y (m)", color=:balance, fmt=:png, colorbar=false)

    ## scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 # convert fractional indices to degree
    py = Array(model.individuals.phytos.sp1.data.y) .* 1 # convert fractional indices to degree
    Plots.scatter!(fl_plot, px, py, ms=5, color = :red, legend=:none)

    ## passive tracer field
    trac1 = Plots.contourf(xC, yC, Array(model.tracers.DOC.data)[3:130,3:130,3]', xlabel="x (m)", ylabel="y (m)", clims=(0.5, 2.4), fmt=:png)

    ## arrange plots side-by-side
    plt = Plots.plot(fl_plot, trac1, size=(800, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" " concentration (mmolC/L)"])

    return plt
end

# The initial condition of the model is below:

plot_model(model)

# The animation is generated using the following code.

begin
    anim = @animate for i in 1:240
       sim_run!(sim)
       plot_model(model)
    end
    gif(anim, "horizontal_2D_example.gif"; fps = 15) # hide
    nothing # hide
end

# ![Simulation animation](horizontal_2D_example.gif)
