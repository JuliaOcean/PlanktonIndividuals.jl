### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 96ebd146-8f33-4e91-a367-97f7ae34a2dc
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()
end

# ╔═╡ 54683b9a-bf00-4956-8894-e263eded3db8
using PlanktonIndividuals, Plots

# ╔═╡ e74d8d84-546b-4ccd-9724-78ed29f587c3
begin
	p=dirname(pathof(PlanktonIndividuals))
	include(joinpath(p,"../examples/helper_functions.jl"))
	nothing
end

# ╔═╡ a1e895f9-fd5a-4908-a760-344b7f8fb304
md"""
# Horizontal 2-Dimensional Example
Here we simulate phytoplankton cells as Lagrangian particles in a horizontal, two-dimensional, flow field.
The domain is periodic in both directions.
"""

# ╔═╡ a6072dc8-2d64-11ec-0a4e-dd442a0f0e63
md"""
## 1. Import packages
"""

# ╔═╡ 88c6027d-bc47-4bc4-9628-c946ac2fd57f
md"""
## 2. Generate Flow Fields
First we generate grid information and the computational architecture (CPU).
"""

# ╔═╡ 10d290fd-5aa7-4b7f-be91-3cc51e4250be
arch = CPU()

# ╔═╡ 3093a948-630d-402c-aeb4-73a552e0cb31
grid = RegularRectilinearGrid(size=(128, 128, 1), spacing=(1, 1, 1))

# ╔═╡ 02a7cda3-4a3a-49ea-a646-e1da7a1e1f14
md"""
Then we use a stream function to generate the flow field which is a double-gyre configuration [as explained here](https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html#Sec7.1)
"""

# ╔═╡ ed5b3467-8bac-4800-824f-2bf26b3954c0
(uvels, vvels, wvels, ϕcenters) = streamfunction_xy();

# ╔═╡ 5abbed9f-a256-4312-a9de-7277379bb33a
md"""
## 3. Model Setup

Next we setup the individual-based model by specifying the architecture, grid, and plankton community.
"""

# ╔═╡ 7735454f-bc9f-4dad-8c40-2e5bba096307
model = PlanktonModel(arch, grid; N_species = 1, 
								  N_individual = 2^7,
								  max_individuals = 2^7*8)

# ╔═╡ 970cba93-d4a5-4579-a436-465085ea6e6b
md"""
Finally we setup the duration of the model simulation and the kind of output we want.
"""

# ╔═╡ e47e4afc-699c-4d6f-957d-26aafb592031
sim = PlanktonSimulation(model, ΔT = 60, iterations = 1,
								vels=(u=uvels, v=vvels, w=wvels),
								ΔT_vel=60*120)

# ╔═╡ ec89f0af-1690-49f5-8d7f-23971d332d7d
md"""
## 4. Model Run

We run the model for 120 time steps (2 hours) and then plot individuals and nutrients in their final state (stored in model).
"""

# ╔═╡ 20a9dffd-b197-44c1-b973-bedea6fd19eb
for i in 1:120
    update!(sim)
end

# ╔═╡ 7f06e8a4-eada-4607-b5dc-8945b8218ca7
md"""
To plot the distribution of individuals as well as nutrient fields we use Plots.jl and create a function that can easily be re-used e.g. to create an animation.
"""

# ╔═╡ f83b3a52-42bb-481e-b405-20f5d4e75f86
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

# ╔═╡ e873d085-a34c-4fd6-a419-5c05444fe0a6
plot_model(model)

# ╔═╡ 9e833367-c0d5-45f0-92a3-0de9042597ff
md"""
Or you can use the following code to generate an animation like below

```
anim = @animate for i in 1:120
   update!(sim)
   plot_model(model)
end
gif(anim, "anim_fps15.gif", fps = 15)
```
![animation](https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/examples/figures/anim_horizontal_2D.gif)
"""

# ╔═╡ Cell order:
# ╟─a1e895f9-fd5a-4908-a760-344b7f8fb304
# ╟─a6072dc8-2d64-11ec-0a4e-dd442a0f0e63
# ╠═96ebd146-8f33-4e91-a367-97f7ae34a2dc
# ╠═54683b9a-bf00-4956-8894-e263eded3db8
# ╟─e74d8d84-546b-4ccd-9724-78ed29f587c3
# ╠═88c6027d-bc47-4bc4-9628-c946ac2fd57f
# ╠═10d290fd-5aa7-4b7f-be91-3cc51e4250be
# ╠═3093a948-630d-402c-aeb4-73a552e0cb31
# ╟─02a7cda3-4a3a-49ea-a646-e1da7a1e1f14
# ╠═ed5b3467-8bac-4800-824f-2bf26b3954c0
# ╟─5abbed9f-a256-4312-a9de-7277379bb33a
# ╠═7735454f-bc9f-4dad-8c40-2e5bba096307
# ╠═970cba93-d4a5-4579-a436-465085ea6e6b
# ╠═e47e4afc-699c-4d6f-957d-26aafb592031
# ╟─ec89f0af-1690-49f5-8d7f-23971d332d7d
# ╠═20a9dffd-b197-44c1-b973-bedea6fd19eb
# ╟─7f06e8a4-eada-4607-b5dc-8945b8218ca7
# ╟─f83b3a52-42bb-481e-b405-20f5d4e75f86
# ╠═e873d085-a34c-4fd6-a419-5c05444fe0a6
# ╟─9e833367-c0d5-45f0-92a3-0de9042597ff
