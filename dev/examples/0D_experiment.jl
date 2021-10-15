### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ b662f1c1-9826-4fe0-9a78-11d1b92386ff
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()
end

# ╔═╡ 1149e16f-752a-4c33-8892-b7c8fa5a2f77
begin
	using PlanktonIndividuals, Plots, JLD2
	using Plots.PlotMeasures
end

# ╔═╡ b795ac78-2d4f-11ec-2cbd-d38135481545
md"""
# Simple 0-Dimensional Example
Here we simulate phytoplankton cells as individuals in a well mixed reactor, like a lab experiment.
"""

# ╔═╡ da41f36d-3747-4edd-8746-4703aed02574
md"""
## 1. Import packages
"""

# ╔═╡ 96fc88b3-f5f4-4882-8757-1ea392610b97
md"""
## 2. Grid Setup
First we generate grid information (one grid box, 256m thick, and 128x128 in width) and the computational architecture (CPU).
"""

# ╔═╡ 28d030b0-8bf7-485d-94f9-5edfffdfd691
arch=CPU()

# ╔═╡ 91079b4b-dd70-4e57-b7ca-4d7be300911a
grid = RegularRectilinearGrid(size=(1,1,1), spacing=(128meters, 128meters, 256meters))

# ╔═╡ c4af7de9-7e1e-44a4-a25c-84cca7c768cf
md"""
## 3. Model Setup
Next we setup the individual-based model by specifying the computational architecture, grid, and plankton community.
"""

# ╔═╡ f7644486-7546-4716-9bb1-58c2220e38eb
model = PlanktonModel(arch, grid; N_species = 1, 
								  N_individual = 2^10, 
								  max_individuals = 2^10*16)

# ╔═╡ e41f726a-3779-47a5-a3ad-ad9a34280f43
md"""
And we setup diagnostics.
"""

# ╔═╡ c199e085-dec2-45e0-88a1-ac409604ad00
diags = PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                                   plankton = (:num, :graz, :mort, :dvid, :PS, :BS),
                                   time_interval = 60seconds)

# ╔═╡ 16295dcc-febe-4350-83e5-e563f51077d1
md"""
Then we setup the duration of the model simulation, a run directory location, and the kind of output we want.
"""

# ╔═╡ 7d99ad3a-79d6-4ede-b25d-e4ac6d03c92e
sim = PlanktonSimulation(model, ΔT = 60seconds, 
								iterations = 1440seconds, 
								diags = diags)

# ╔═╡ d7852aae-c55b-490b-8b69-2594d4354810
md"""
Finally we setup the output writer.
"""

# ╔═╡ 0cc1d551-d419-4573-907b-dadfdb39e2da
sim.output_writer = PlanktonOutputWriter(save_diags = true)

# ╔═╡ 0fcd4f42-1ac4-4a0e-86f6-16e6454246d2
md"""
## 4. Model Run
"""

# ╔═╡ 632bfd5d-2879-4669-88f7-779c5ff6aa78
update!(sim)

# ╔═╡ 3ddc5181-28fe-41d9-9947-0f059d99709e
md"""
## 5. Results Vizualization
"""

# ╔═╡ c41e9408-89e6-47ad-93cd-7d65bb173972
md"""
Open the output file
"""

# ╔═╡ dfad80f2-375e-4264-9ea4-5ad69ef9a766
file = jldopen(sim.output_writer.diags_file, "r")

# ╔═╡ 21e59e96-2a01-41a6-b1dd-234048501f23
md"""
Extract a vector of iterations
"""

# ╔═╡ e9224442-e4ec-4475-a9af-aadd42a8e8e5
iterations = parse.(Int, keys(file["timeseries/t"]))

# ╔═╡ f436d0fa-1427-4da7-a28e-ec596c818546
md"""
Read results into arrays and close file.
"""

# ╔═╡ 4d09f435-976b-4b86-b0b7-7931bf4bd99d
begin
	num  = zeros(1440) 
	dvid = zeros(1440)
	mort = zeros(1440)
	PS   = zeros(1440)
	for (i, iter) in enumerate(iterations)
		num[i]  = file["timeseries/num/$iter"][1,1,1]
		dvid[i] = file["timeseries/dvid/$iter"][1,1,1]
		mort[i] = file["timeseries/mort/$iter"][1,1,1]
		PS[i]   = file["timeseries/PS/$iter"][1,1,1]
	end
	close(file)
end

# ╔═╡ 49c7d833-6fe1-4619-ace3-e2f6aaefac32
md"""
Now we plot the results
"""

# ╔═╡ d642e1ee-7e66-46d6-9f92-8b29279e7ad9
begin
	p1 = plot(collect(1:60:60*1440) ./ 3600, num, title = "population", legend=:none, fmt=:png, bottom_margin = 5mm)
	p2 = plot(collect(1:60:60*1440) ./ 3600, dvid ./ num .* 60, title = "division rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
	p3 = plot(collect(1:60:60*1440) ./ 3600, mort ./ num .* 60, title = "mortarlity rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
	p4 = plot(collect(1:60:60*1440) ./ 3600, PS ./ num .* 12 .* 1e12 .* 3600, title = "photosynthesis rate (fg C/cell/hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
	nothing
end

# ╔═╡ f3ac2cf9-5b6c-453f-ac2f-f082277e0cbb
plot(p1, p2, p3, p4, layout = (4,1), size=(600,600), titlefont = (12))

# ╔═╡ Cell order:
# ╟─b795ac78-2d4f-11ec-2cbd-d38135481545
# ╠═b662f1c1-9826-4fe0-9a78-11d1b92386ff
# ╟─da41f36d-3747-4edd-8746-4703aed02574
# ╠═1149e16f-752a-4c33-8892-b7c8fa5a2f77
# ╟─96fc88b3-f5f4-4882-8757-1ea392610b97
# ╠═28d030b0-8bf7-485d-94f9-5edfffdfd691
# ╠═91079b4b-dd70-4e57-b7ca-4d7be300911a
# ╟─c4af7de9-7e1e-44a4-a25c-84cca7c768cf
# ╠═f7644486-7546-4716-9bb1-58c2220e38eb
# ╟─e41f726a-3779-47a5-a3ad-ad9a34280f43
# ╠═c199e085-dec2-45e0-88a1-ac409604ad00
# ╟─16295dcc-febe-4350-83e5-e563f51077d1
# ╠═7d99ad3a-79d6-4ede-b25d-e4ac6d03c92e
# ╟─d7852aae-c55b-490b-8b69-2594d4354810
# ╠═0cc1d551-d419-4573-907b-dadfdb39e2da
# ╟─0fcd4f42-1ac4-4a0e-86f6-16e6454246d2
# ╠═632bfd5d-2879-4669-88f7-779c5ff6aa78
# ╟─3ddc5181-28fe-41d9-9947-0f059d99709e
# ╟─c41e9408-89e6-47ad-93cd-7d65bb173972
# ╠═dfad80f2-375e-4264-9ea4-5ad69ef9a766
# ╟─21e59e96-2a01-41a6-b1dd-234048501f23
# ╠═e9224442-e4ec-4475-a9af-aadd42a8e8e5
# ╟─f436d0fa-1427-4da7-a28e-ec596c818546
# ╠═4d09f435-976b-4b86-b0b7-7931bf4bd99d
# ╟─49c7d833-6fe1-4619-ace3-e2f6aaefac32
# ╠═d642e1ee-7e66-46d6-9f92-8b29279e7ad9
# ╠═f3ac2cf9-5b6c-453f-ac2f-f082277e0cbb
