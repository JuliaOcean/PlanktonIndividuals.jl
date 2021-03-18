# # Vertical 2-Dimensional Example
#
# Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field.
# The domain is periodic in x direction but bounded in z direction

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we'll generate grid information
grid = gen_Grid(size=(128, 1, 128), spacing=(1, 1, 1))

# Then we use a stream function to generate the flow fields

f(x, y, z) = cos(x*2π/128)*cos(z*2π/128) #stream function
u(x, y, z) = -sin(z*2π/128)*2π/128
w(x, y, z) = sin(x*2π/128)*2π/128
scal = 2e-2
ϕ  = [f(x, y, z) for x in grid.xC[3:130], y in grid.yC[3:3], z in grid.zC[3:130]] .* scal

uF = [u(x, y, z) for x in grid.xF[3:131], y in grid.yF[3:3], z in grid.zF[3:131]] .* scal # array size is (129,1,129)
wF = [w(x, y, z) for x in grid.xF[3:131], y in grid.yF[3:3], z in grid.zF[3:131]] .* scal # array size is (129,1,129)
vF = zeros(128,1,128)

uvels = fill(uF[1:128,1:1,1:128], 2)
vvels = fill(vF, 2)
wvels = fill(wF[1:128,1:1,1:129], 2)
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

sim = PI_simulation(model, ΔT = 60, nΔT = 1, diag_freq = 3600, 
                    vels=(u=uvels, v=vvels, w=wvels), 
                    vel_reuse = true)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Run the Model
#
# Finally, we run the model for 60 time steps (1 hour) and plot the individuals and DOC field.

for i in 1:60
    update!(sim)
end

## Coordinate arrays for plotting
xC, zC = collect(model.grid.xC)[3:130], collect(model.grid.zC)[3:130]

## heatmap of the flow field
fl_plot = heatmap(xC, zC, ϕ', xlabel="x (m)", ylabel="z (m)", color=:balance, fmt=:png)

## a scatter plot embeded in the flow fields
px = Array(model.individuals.phytos.sp1.data.x)
pz = Array(model.individuals.phytos.sp1.data.z)
Plots.scatter!(fl_plot, px, pz, ms=5, color = :red, legend=:none)


#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# Or you can use the following code to generate an animation
#
# ```
# anim = @animate for i in 1:100
#     update!(sim)

#   ## Coordinate arrays for plotting
#     xC, zC = collect(model.grid.xC)[3:130], collect(model.grid.zC)[3:130]

#     ## heatmap of the flow field
#     fl_plot = heatmap(xC, zC, ϕ[:,1,:], xlabel="x (m)", ylabel="y (m)", color=:balance, clims=(-5e-2, 5e-2))

#     ## a scatter plot embeded in the flow fields
#     px = Array(model.individuals.phytos.sp1.data.x)
#     pz = Array(model.individuals.phytos.sp1.data.z)
#     Plots.scatter!(fl_plot, px, pz, ms=5, color = :red, legend=:none)

#     ## DOC field
#     trac1 = Plots.heatmap(xC, zC, Array(model.nutrients.DOC.data)[3:130,3,3:130]', clims=(0.10,1.05), xlabel="x (m)", ylabel="z (m)")

#     ## Arrange the plots side-by-side.
#     plot(fl_plot, trac1, size=(1200, 400),
#          title=[lpad(i÷1440,2,"0")*"day "*lpad(i÷60-24*(i÷1440),2,"0")*"hour" "DOC (mmolC/L)"])
# end
# gif(anim, "anim.gif", fps = 15)
# ```