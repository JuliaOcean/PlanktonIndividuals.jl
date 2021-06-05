# # Vertical 2-Dimensional Example
#
# Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field, with 
# one horizontal direction (x) and one vertical one (z), like in an ocean transect.
#
# Here the domain is periodic in the x direction while it is bounded in the z direction.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots

p=dirname(pathof(PlanktonIndividuals))
include(joinpath(p,"../examples/helper_functions.jl"))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we generate grid information (128 by 128 grid boxes, 1m thick, and 1m wide) and the computational architecture (CPU).

arch = CPU()

grid = RegularRectilinearGrid(size=(128, 1, 128), spacing=(1, 1, 1))

# Then we use a stream function (see helper_functions.jl) to generate a simple flow field (displayed below)
# in a 2D vertical plane.

(uvels, vvels, wvels, ϕcenters) = streamfunction_xz();

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid, and plankton community.

model = PlanktonModel(arch, grid; N_species = 1, N_individual = 2^7, max_individuals = 2^7*8)

# Finally we setup the duration of the model simulation and the kind of output we want.

sim = PlanktonSimulation(model, ΔT = 60, nΔT = 1, vels=(u=uvels, v=vvels, w=wvels), vel_reuse = true)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Model Run
#
# We run the model for 120 time steps (2 hours) and then plot individuals and nutrients in their final state (stored in model).

for i in 1:120
    update!(sim)
end

# To plot the distribution of individuals as well as nutrient fields we use Plots.jl and 
# create a function that can easily be re-used e.g. to create an animation.

function plot_model(model::PlanktonModel)
    ## Coordinate arrays for plotting
    xC, zC = collect(model.grid.xC)[3:130], collect(model.grid.zC)[3:130]

    ## contour of the flow field
    fl_plot = Plots.contourf(xC, reverse(zC), rotl90(ϕcenters), xlabel="x (m)", ylabel="z (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 # convert fractional indices to degree
    pz = Array(model.individuals.phytos.sp1.data.z) .* -1# convert fractional indices to degree
    Plots.scatter!(fl_plot, px, pz, ms=5, color = :red, legend=:none)

    ## DOC field
    trac1 = Plots.contourf(xC, reverse(zC), rotl90(Array(model.nutrients.DOC.data)[3:130,3,3:130]), xlabel="x (m)", ylabel="z (m)", clims=(0.5, 1.1), fmt=:png)

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
# ![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_vertical_2D.gif)