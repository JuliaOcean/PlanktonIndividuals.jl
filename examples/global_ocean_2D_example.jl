# #  Global Ocean Example
#
# Here we simulate phytoplankton cells as **passive** Lagrangian particles and nutrient fields as **passive** tracers in the global ocean.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots, IndividualDisplacements, MeshArrays, OceanStateEstimation
pyplot()

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Generate Flow Fields
#
# First we'll generate grid information
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
IndividualDisplacements.get_occa_velocity_if_needed()
𝑃,𝐷,Γ=OCCA_FlowFields()

grid = RegularLatLonGrid(size=(360, 160, 1), lat=(-80,80), lon=(-180,180), z=(0,-10))

# Next, we generate a mask from  the land shape information 
landshape = findall(x -> x == 0.0 ,Γ.Depth[1])
mask = ones(360,160)
mask[landshape] .= 0.0
mask = reshape(mask,360,160,1)
nothing

# Then we re-format velocity fields so they can be loaded in to PlanktonSimulation
uu=reshape(𝑃.u0[1][2:end-1,2:end-1,1],(360,160,1)) .* grid.dxF[3:end-2, 3:end-2] .* mask
vv=reshape(𝑃.v0[1][2:end-1,2:end-1,1],(360,160,1)) .* grid.dyF[3:end-2, 3:end-2] .* mask
vv = hcat(vv,zeros(360)) # for bounded boundary condition
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

#
# In order to change individuals and tracer to **passive** mode, we need to update some parameters.
#
phyt_parameter = Dict("PCmax"   => [0.0], # maximum photosynthesis rate
                      "VNH4max" => [0.0], # maximum ammonia uptake rate
                      "VNO3max" => [0.0], # maximum nitrate uptake rate
                      "VPO4max" => [0.0], # maximum phosphate uptake rate
                      "respir_a"=> [0.0], # respiration rate
                      "k_mtb"   => [0.0], # biosynthesis rate
                      "dvid_P"  => [0.0], # probability of cell division
                      "mort_P"  => [0.0]  # probability of cell natural death
                     )
#

bgc_parameter = Dict("kDOC"     => 0.0,         # Remineralization rate for DOC, turn over time: a month (per second)
                     "Nit"      => 0.0,         # Nitrification rate for NH4
                     "kDON"     => 0.0,         # Remineralization rate for DON, turn over time: a month (per second)
                     "kDOP"     => 0.0,         # Remineralization rate for DON, turn over time: a month (per second)
                     "kPOC"     => 0.0,         # Remineralization rate for POC, turn over time: a month (per second)
                     "kPON"     => 0.0,         # Remineralization rate for PON, turn over time: a month (per second)
                     "kPOP"     => 0.0,         # Remineralization rate for PON, turn over time: a month (per second)
                    )
#
# Then we update new parameter values in the model
model = PlanktonModel(CPU(), grid;
                      N_species = 1,
                      N_individual = 360,
                      max_individuals = 360*8,
                      bgc_params = update_bgc_params(bgc_parameter),
                      phyt_params = update_phyt_params(phyt_parameter), 
                      mask = mask)

# We also need to setup a runtime simulation to run the model.
# The simulation includes time step, number of time steps, flow fields that
# will be used etc.

sim = PlanktonSimulation(model, ΔT = 3600, iterations = 1, vels=(u=uvels, v=vvels, w=wvels), ΔT_vel=3600*24)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Run the Model
#
# Finaly, we run the model and plot the distribution of individuals as well as nutrient fields
# We use Plots.jl to plot individuals and nutrient fields.
#
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
    trac1 = Plots.heatmap(xC, yC, Array(model.nutrients.DOC.data)[3:end-2,3:end-2,3]', xlabel="x (m)", ylabel="y (m)", clims=(0.5, 1.1), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(1200, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "DOC (mmolC/L)"])

    return plt
end
#
# We run the model for 24 time steps (1 hour per time step) and plot the individuals and DOC field.
for i in 1:24
    update!(sim)
end

#
# We plot the current state of the model
u_plot = uu[:,:,1]
u_plot[landshape] .= NaN
plot_model(model, u_plot)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# Or you can use the following code to generate an animation like below.
# Please note that the following simulation is run for a year.
#
# ```
# anim = @animate for i in 1:120
#   update!(sim)
#   plot_model(model)
# end
# gif(anim, "anim_fps15.gif", fps = 15)
# ```
# ![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_global.gif)