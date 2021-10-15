### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 5bd797cb-5042-4cc5-b988-48d728033156
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()
end

# ╔═╡ 209bc3b1-e354-49da-9f49-a76f6a82c445
using PlanktonIndividuals, Plots, JLD2

# ╔═╡ 6108b4f8-2d69-11ec-20c7-4d58443806cf
md"""
#  Global Ocean Example

Here we simulate phytoplankton cells as **passive** Lagrangian particles and nutrient fields as **passive** tracers in the global ocean.
"""

# ╔═╡ 88855dc0-ccc7-49b9-b066-4c9059018c36
md"""
## 1. Import packages
"""

# ╔═╡ b523e344-f0d4-424e-82a2-f571e128d9ab
md"""
## 2. Generate Flow Fields

First we'll generate grid information
"""

# ╔═╡ 707f61f5-8b6a-4c2b-8b5b-3aedd62e0f9c
grid = RegularLatLonGrid(size=(360, 160, 1), lat=(-80,80), lon=(-180,180), z=(0,-10))

# ╔═╡ 4871411e-1d61-46dc-81ff-8143c14fb924
md"""
Next, we generate a mask from  the land shape information.
"""

# ╔═╡ 2451216d-cba2-49fe-acc6-91673dae48ad
md"""
Then we re-format velocity fields so they can be loaded in to PlanktonSimulation.
"""

# ╔═╡ 7582dc85-d30a-4996-b348-526e44d9b780
begin
	file = jldopen(PlanktonIndividuals.global_vels,"r")
	mask = file["mask"]
	uvels = file["u"]
	vvels = file["v"]
	wvels = file["w"]
	close(file)
end

# ╔═╡ d88cd912-326b-4d21-ac18-74e3ff7edcb9
md"""
## 3. Model Setup

Next we setup the individual-based model by specifying the architecture, grid, number of individuals, parameters, and nutrient initial conditions.

In order to change individuals and tracer to **passive** mode, we need to update some parameters.
"""

# ╔═╡ edfb412d-d002-43fd-9575-0a86b16be12e
phyt_parameter = Dict("PCmax"   => [0.0], # maximum photosynthesis rate
                      "respir_a"=> [0.0], # respiration rate
                      "dvid_P"  => [0.0], # probability of cell division
                      "mort_P"  => [0.0]  # probability of cell natural death
                     )


# ╔═╡ 2fed20a6-075d-4247-99ec-b7ec9797280b
bgc_parameter = Dict("kDOC" => 0.0, # Remineralization rate for DOC
                     "Nit"  => 0.0, # Nitrification rate for NH4
                     "kDON" => 0.0, # Remineralization rate for DON
                     "kDOP" => 0.0, # Remineralization rate for DON
                     "kPOC" => 0.0, # Remineralization rate for POC
                     "kPON" => 0.0, # Remineralization rate for PON
                     "kPOP" => 0.0, # Remineralization rate for PON
                    )

# ╔═╡ d085c56b-b53a-404d-a7cd-85ab7b9dca54
md"""
Then we update new parameter values in the model.
"""

# ╔═╡ 1fb28de1-536c-4fe1-bcfb-99871e7f25fd
model = PlanktonModel(CPU(), grid;
					  mode = CarbonMode(),
                      N_species = 1,
                      N_individual = 360,
                      max_individuals = 360*8,
                      bgc_params = update_bgc_params(bgc_parameter),
                      phyt_params = update_phyt_params(phyt_parameter, CarbonMode()), 
                      mask = mask)

# ╔═╡ b1b650bf-c997-4bba-a033-15e5897e1479
md"""
We also need to setup a runtime simulation to run the model.
The simulation includes time step, number of time steps, flow fields that
will be used etc.
"""

# ╔═╡ 6f395d24-20ca-480a-a0db-95f7e1f37a0b
sim = PlanktonSimulation(model, ΔT = 3600,
								iterations = 1,
								vels=(u=uvels, v=vvels, w=wvels),
								ΔT_vel=3600*24)

# ╔═╡ 0513073a-196b-4bea-8f15-baa6f8b20c7f
md"""
## 4. Run the Model

Finaly, we run the model and plot the distribution of individuals as well as nutrient fields.

We use Plots.jl to plot individuals and nutrient fields.
"""

# ╔═╡ 01d4a913-2f06-49f3-b241-22e07acd6b3f
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

# ╔═╡ 01ddf7ae-7609-4dcf-99dc-0622e3960a56
md"""
We run the model for 24 time steps (1 hour per time step) and plot the individuals and DOC field.
"""

# ╔═╡ da402bd4-797b-4dc1-8c9b-ac930f760754
for i in 1:24
    update!(sim)
end

# ╔═╡ 03152520-d23e-4d5c-8cc1-5445e296799c
md"""
We plot the current state of the model.
"""

# ╔═╡ 2dfe1e64-bd04-457e-8a85-bf60cf0c724e
begin
	u_plot = uvels[:,:,1,1]
	u_plot[findall(x -> x == 0.0, mask)] .= NaN
	plot_model(model, u_plot)
end

# ╔═╡ f5c2f5ab-e875-4904-98a2-162ec52d700e
md"""
Or you can use the following code to generate an animation like below.

Please note that the following simulation is run for a year.

```
anim = @animate for i in 1:120
  update!(sim)
  plot_model(model)
end
gif(anim, "anim_fps15.gif", fps = 15)
```
![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_global.gif)
"""

# ╔═╡ Cell order:
# ╟─6108b4f8-2d69-11ec-20c7-4d58443806cf
# ╟─88855dc0-ccc7-49b9-b066-4c9059018c36
# ╠═5bd797cb-5042-4cc5-b988-48d728033156
# ╠═209bc3b1-e354-49da-9f49-a76f6a82c445
# ╟─b523e344-f0d4-424e-82a2-f571e128d9ab
# ╠═707f61f5-8b6a-4c2b-8b5b-3aedd62e0f9c
# ╟─4871411e-1d61-46dc-81ff-8143c14fb924
# ╟─2451216d-cba2-49fe-acc6-91673dae48ad
# ╠═7582dc85-d30a-4996-b348-526e44d9b780
# ╟─d88cd912-326b-4d21-ac18-74e3ff7edcb9
# ╠═edfb412d-d002-43fd-9575-0a86b16be12e
# ╠═2fed20a6-075d-4247-99ec-b7ec9797280b
# ╟─d085c56b-b53a-404d-a7cd-85ab7b9dca54
# ╠═1fb28de1-536c-4fe1-bcfb-99871e7f25fd
# ╟─b1b650bf-c997-4bba-a033-15e5897e1479
# ╠═6f395d24-20ca-480a-a0db-95f7e1f37a0b
# ╟─0513073a-196b-4bea-8f15-baa6f8b20c7f
# ╟─01d4a913-2f06-49f3-b241-22e07acd6b3f
# ╟─01ddf7ae-7609-4dcf-99dc-0622e3960a56
# ╠═da402bd4-797b-4dc1-8c9b-ac930f760754
# ╟─03152520-d23e-4d5c-8cc1-5445e296799c
# ╠═2dfe1e64-bd04-457e-8a85-bf60cf0c724e
# ╟─f5c2f5ab-e875-4904-98a2-162ec52d700e
