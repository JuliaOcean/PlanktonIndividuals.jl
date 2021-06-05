# # Horizontal 2-Dimensional Example
#
# Here we simulate phytoplankton cells as Lagrangian particles in a horizontal, two-dimensional, flow field.
# The domain is periodic in both directions.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots

p=dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"../examples/helper_functions.jl"))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we'll generate grid information
grid = RegularRectilinearGrid(size=(128, 128, 1), spacing=(1, 1, 1))

# Then we use a stream function to generate the flow field which is a double-gyre configuration
# [as explained here](https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html#Sec7.1)

(uvels, vvels, wvels, ϕcenters) = streamfunction_xy();

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Set Up The Model
#
# Next we setup the individual-based model by specifying the architecture, grid,
# number of individuals, parameters, and nutrient initial conditions.

model = PlanktonModel(CPU(), grid; N_species = 1, N_individual = 2^7, max_individuals = 2^7*8)

# We also need to setup a runtime simulation to run the model.
# The simulation includes time step, number of time steps, flow fields that
# will be used etc.

sim = PlanktonSimulation(model, ΔT = 60, nΔT = 1, vels=(u=uvels, v=vvels, w=wvels), vel_reuse = true)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Run the Model
#
# We run the model for 120 time steps (2 hour) and plot the individuals and DOC field afterwards.
for i in 1:120
    update!(sim)
end

# To plot the distribution of individuals as well as nutrient fields we use Plots.jl and 
# create a function that can easily be re-used e.g. to create an animation.

function plot_model(model::PlanktonModel)
    ## Coordinate arrays for plotting
    xC, yC = collect(model.grid.xC)[3:130], collect(model.grid.yC)[3:130]

    ## heatmap of the flow field
    fl_plot = Plots.contourf(xC, yC, ϕcenters', xlabel="x (m)", ylabel="y (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 # convert fractional indices to degree
    py = Array(model.individuals.phytos.sp1.data.y) .* 1 # convert fractional indices to degree
    Plots.scatter!(fl_plot, px, py, ms=5, color = :red, legend=:none)

    ## DOC field
    trac1 = Plots.contourf(xC, yC, Array(model.nutrients.DOC.data)[3:130,3:130,3]', xlabel="x (m)", ylabel="y (m)", clims=(0.5, 1.1), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(800, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "DOC (mmolC/L)"])

    return plt
end

plot_model(model)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# Or you can use the following code to generate an animation like below
#
# ```
# anim = @animate for i in 1:120
#    update!(sim)
#    plot_model(model)
# end
# gif(anim, "anim_fps15.gif", fps = 15)
# ```
# ![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_horizontal_2D.gif)