# # Simple Lab Example
#
# Here we simulate phytoplankton cells as individuals in a well mixed reactor, like a lab experiment.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import packages
#
using PlanktonIndividuals, Plots, JLD2
using Plots.PlotMeasures

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Grid Setup
#
# First we generate grid information (one grid box, 256m thick, and 128x128 in width) and the computational architecture (CPU).
grid = RegularRectilinearGrid(size=(1,1,1), spacing=(128, 128, 256))

arch=CPU()

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Model Setup
#
# Next we setup the individual-based model by specifying the computational architecture, grid, plankton community, and diagnostics.
#
model = PI_Model(arch, grid; 
                 individual_size = (Nsp = 1, N = 2^10, cap = 16),
                 diag_nprocs = (:num, :graz, :mort, :dvid, :PS, :BS, :Cq, :chl, :Bm))

# Finally we setup the duration of the model simulation, a run directory location, and the kind of output we want.
#
ntimesteps = 60 * 24 * 1 # 1 days with time step of 60s 
res_dir = PrepRunDir()
sim = PI_simulation(model, ΔT = 60, nΔT = ntimesteps, diag_freq = 1, 
                    res_dir = res_dir, 
                    save_diags = true,
                    save_individuals = false)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Model Run
#
update!(sim)


#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 5. Results Vizualization
#
# Empty arrays to store results:
#
num  = zeros(1440) 
dvid = zeros(1440)
mort = zeros(1440)
PS   = zeros(1440)
nothing
#
# Read results from output file:
#
path = joinpath(res_dir,"diags.jld2")
file = jldopen(path, "r")
for i in 1:1440
    num[i]  = file[joinpath(lpad(i*60,10,"0"),"sp1","num")][3,3,3]
    dvid[i] = file[joinpath(lpad(i*60,10,"0"),"sp1","dvid")][3,3,3]
    mort[i] = file[joinpath(lpad(i*60,10,"0"),"sp1","mort")][3,3,3]
    PS[i]   = file[joinpath(lpad(i*60,10,"0"),"sp1","PS")][3,3,3]
end
close(file)
#
# Plot results:
#
p1 = plot(collect(1:60:60*1440) ./ 3600, num, title = "population", legend=:none, fmt=:png, bottom_margin = 5mm)
p2 = plot(collect(1:60:60*1440) ./ 3600, dvid ./ num .* 60, title = "division rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
p3 = plot(collect(1:60:60*1440) ./ 3600, mort ./ num .* 60, title = "mortarlity rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
p4 = plot(collect(1:60:60*1440) ./ 3600, PS ./ num .* 12 .* 1e12 .* 3600, title = "photosynthesis rate (fg C/cell/hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
plot(p1, p2, p3, p4, layout = (4,1), size=(600,600), titlefont = (12))