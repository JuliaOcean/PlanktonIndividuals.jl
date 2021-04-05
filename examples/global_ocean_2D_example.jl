# # Horizontal 2-Dimensional Example
#
# Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field.
# The domain is periodic in both directions horizontally.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots, IndividualDisplacements, MeshArrays, OceanStateEstimation

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we'll generate grid information
grid = RegularRectilinearGrid(size=(360, 160, 1), spacing=(1, 1, 1))

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
IndividualDisplacements.get_occa_velocity_if_needed()
𝑃,𝐷,Γ=OCCA_FlowFields(backward_in_time=false)

# Then we use a stream function to generate the flow field which is a double-gyre configuration
uu=reshape(𝑃.u0[1][2:end-1,2:end-1,1],(360,160,1))
vv=reshape(𝑃.v0[1][2:end-1,2:end-1,1],(360,160,1))
ww=zeros(360,160,2)

uvels = fill(uu, 2)
vvels = fill(vv, 2)
wvels = fill(ww, 2)
uvels = cat(uvels..., dims=4)
vvels = cat(vvels..., dims=4)
wvels = cat(wvels..., dims=4)

nothing

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid,
# number of individuals, parameters, and nutrient initial conditions.

model = PI_Model(CPU(), grid; individual_size = (Nsp = 1, N = 360, cap = 8),
                 nut_source = [1.0, 0.02, 0.05, 0.01, 1.0, 0.1, 0.02, 0.2, 0.02, 0.001])

# We also need to setup a runtime simulation to run the model.
# The simulation includes time step, number of time steps, flow fields that
# will be used etc.

sim = PI_simulation(model, ΔT = 3600, nΔT = 1, diag_freq = 86400, 
                    vels=(u=uvels, v=vvels, w=wvels), 
                    vel_reuse = true)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Run the Model
#
# Finaly, we run the model and plot the distribution of individuals as well as nutrient fields
# We use Plots.jl to plot individuals and nutrient fields.
#
function plot(model::PI_Model)
    ## Coordinate arrays for plotting
    xC, yC = collect(model.grid.xC)[3:end-2], collect(model.grid.yC)[3:end-2]

    ## heatmap of the flow field
    fl_plot = Plots.contourf(xC, yC, uu[:,:,1]', xlabel="x (m)", ylabel="y (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x)
    py = Array(model.individuals.phytos.sp1.data.y)
    Plots.scatter!(fl_plot, px, py, ms=5, color = :red, legend=:none)

    ## DOC field
    trac1 = Plots.contourf(xC, yC, Array(model.nutrients.DOC.data)[3:end-2,3:end-2,3]', xlabel="x (m)", ylabel="y (m)", clims=(0.5, 1.1), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(800, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "DOC (mmolC/L)"])

    return plt
end
#
# We run the model for 24 time steps (1 hour) and plot the individuals and DOC field.
for i in 1:24
    println("$i")
    update!(sim)
end

plot(model)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# Or you can use the following code to generate an animation like below
#
# ```
# anim = @animate for i in 1:120
#   update!(sim)
#   plot(model)
# end
# gif(anim, "anim_fps15.gif", fps = 15)
# ```
# ![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_horizontal_2D.gif)