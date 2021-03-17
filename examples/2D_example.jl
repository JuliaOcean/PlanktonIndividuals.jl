# # Horizontal 2-Dimensional Example
#
# Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field.
# The domain is periodic in both directions Horizontally.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we'll generate grid information
grid = gen_Grid(size=(32, 32, 1), spacing=(1, 1, 1))

# Then we use a stream function to generate the flow fields

f(x, y, z) = sin(x/8) + cos(y/8) #stream function
u(x, y, z) = -sin(y/8)/8
v(x, y, z) = -cos(x/8)/8
scal = 2e-2
ϕ  = [f(x, y, z) for x in grid.xC[3:34], y in grid.yC[3:34], z in grid.zC[3]] .* scal
uC = [u(x, y, z) for x in grid.xC[3:34], y in grid.yC[3:34], z in grid.zC[3]] .* scal
vC = [v(x, y, z) for x in grid.xC[3:34], y in grid.yC[3:34], z in grid.zC[3]] .* scal
wC = zeros(32,32,1)

uF=0.5*(circshift(uC, (1,0))+uC)
vF=0.5*(circshift(vC, (0,1))+vC)
wF=0.5*(circshift(wC, (0,0))+wC)

uvels = fill(uF, 11)
vvels = fill(vF, 11)
wvels = fill(wF, 11)
uvels = cat(uvels..., dims=4)
vvels = cat(vvels..., dims=4)
wvels = cat(wvels..., dims=4)
nothing

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the architecture, grid,
# number of individuals, parameters, and nutrient initial conditions.

model = PI_Model(CPUs(), grid; individual_size = (Nsp = 1, N = 2^7, cap = 8),
                 nut_source = [1.0, 0.02, 0.05, 0.01, 1.0, 0.1, 0.02, 0.2, 0.02, 0.001])

# We also need to setup a runtime simulation to run the model.
# The simulation includes time step, number of time steps, flow fields that
# will be used etc.

sim = PI_simulation(model, ΔT = 60, nΔT = 10, diag_freq = 3600, 
                    vels=(u=uvels, v=vvels, w=wvels), 
                    vel_reuse = true)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Run the Model
#
# Finally, we run the model and plot the individuals and DOC field.

update!(sim)

## Coordinate arrays for plotting
xC, yC = collect(model.grid.xC)[3:34], collect(model.grid.yC)[3:34]

## heatmap of the flow field
fl_plot = heatmap(xC, yC, ϕ', xlabel="x (m)", ylabel="y (m)", color=:balance, xlims=(0,32), ylims=(0,32), clims=(-5e-2, 5e-2))

## a scatter plot embeded in the flow fields
px = Array(model.individuals.phytos.sp1.data.x)
py = Array(model.individuals.phytos.sp1.data.y)
Plots.scatter!(fl_plot, px, py, ms=5, color = :red, legend=:none)

## DOC field
trac1 = Plots.heatmap(xC, yC, Array(model.nutrients.DOC.data)[3:34,3:34,3]', clims=(0.10,1.05), xlabel="x (m)", ylabel="y (m)")

## Arrange the plots side-by-side.
plt = plot(fl_plot, trac1, size=(800, 400), title=["Individuals" "DOC (mmolC/L)"])
# display(plt)


#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# Or you can use the following code to generate an animation
#
# ```
# sim.nΔT = 1 # stop every time step to plot
# anim = @animate for i in 1:100
#     update!(sim)
#   
#     ## Coordinate arrays for plotting
#     xC, yC = collect(model.grid.xC)[3:34], collect(model.grid.yC)[3:34]
#
#     ## heatmap of the flow field
#     fl_plot = heatmap(xC, yC, ϕ', xlabel="x (m)", ylabel="y (m)", color=:balance, xlims=(0,32), ylims=(0,32), clims=(-5e-2, 5e-2))
#
#     ## a scatter plot embeded in the flow fields
#     px = Array(model.individuals.phytos.sp1.data.x)
#     py = Array(model.individuals.phytos.sp1.data.y)
#     Plots.scatter!(fl_plot, px, py, ms=5, color = :red, legend=:none)
#
#     ## DOC field
#     trac1 = Plots.heatmap(xC, yC, Array(model.nutrients.DOC.data)[3:34,3:34,3]', clims=(0.10,1.05), xlabel="x (m)", ylabel="y (m)")
#
#     ## Arrange the plots side-by-side.
#     plot(fl_plot, trac1, size=(1200, 400),
#          title=[lpad(i÷1440,2,"0")*"day "*lpad(i÷60-24*(i÷1440),2,"0")*"hour" "DOC (mmolC/L)"])
# end
#
# gif(anim, "anim.gif", fps = 15)
# ```