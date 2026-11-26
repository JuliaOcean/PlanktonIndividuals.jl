# # Vertical 2-Dimensional Example
# Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field, with one horizontal direction (x) and one vertical one (z), like in an ocean transect.
#
# Here the domain is periodic in the x direction while it is bounded in the z direction.

# ## 1. Import packages

using PlanktonKernels.Units
using PlanktonIndividuals, Plots, PlanktonKernels

include(joinpath(pkgdir(PlanktonIndividuals), "examples", "helper_functions.jl"))

# ## 2. Generate Flow Fields
# First we generate grid information (128 by 128 grid boxes, 1m thick, and 1m wide) and the computational architecture (CPU).

arch = CPU()

grid = RectilinearGrid(size=(128, 1, 128), x = (0.0, 128meters), y = (0.0,1meter), z = (0.0, -128meter))

# Then we use a stream function (see helper_functions.jl) to generate a simple flow field (displayed below) in a 2D vertical plane.

(uvels, vvels, wvels, ϕcenters) = streamfunction_xz();

# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid, and plankton community.

model = PlanktonModel(arch, grid; phyto = phyto_setup(nothing, [2^7], 1),
                                  options = model_options(max_individuals = 2^7*8))

# We modify the DOC tracer field to a two-layer setup to highlight advection.
#
# !!! note
#     DOC is used here as a generic tracer; not to represent dissolved organic carbon per se.

initialized_model = begin
    model.tracers.DOC.data[:,:,1:66] .= 1.5
    model.tracers.DOC.data[:,:,67:end] .= 2.0
    model
end

# Finally we setup the duration of the model simulation and the kind of output we want.

sim = PlanktonSimulation(initialized_model; ΔT = 60.0, iterations = 1,
                                vels=(u=uvels, v=vvels, w=wvels),
                                ΔT_vel=60.0*120)

# ## 4. Model Run
#
# We run the model for 240 time steps (4 hours) and generate an animation.

# To plot the distribution of individuals as well as nutrient fields we use Plots.jl and create a function that can easily be re-used e.g. to create an animation.

function plot_model(model::PlanktonModel)
    ## coordinate arrays for plotting
    xC, zC = collect(model.grid.xC)[3:130], collect(model.grid.zC)[3:130]

    ## contour of the flow field
    fl_plot = Plots.contourf(xC, reverse(zC), rotl90(ϕcenters), xlabel="x (m)", ylabel="z (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 # convert fractional indices to degree
    pz = Array(model.individuals.phytos.sp1.data.z) .* -1# convert fractional indices to degree
    Plots.scatter!(fl_plot, px, pz, ms=5, color = :red, legend=:none)

    ## tracer field
    trac1 = Plots.contourf(xC, reverse(zC), rotl90(Array(model.tracers.DOC.data)[3:130,3,3:130]), xlabel="x (m)", ylabel="z (m)", clims=(0.0, 2.5), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(800, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "DOC (mmolC/L)"])

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
    gif(anim, "vertical_2D_example.gif"; fps = 15) # hide
    nothing # hide
end

# ![Simulation animation](vertical_2D_example.gif)
