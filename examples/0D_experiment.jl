# # Simple 0-Dimensional Example
# Here we simulate phytoplankton cells as individuals in a well mixed reactor, like a lab experiment.

# ## 1. Import packages

using PlanktonKernels.Units
using PlanktonIndividuals, Plots, JLD2, PlanktonKernels
using Plots.PlotMeasures

# ## 2. Grid Setup
# First we generate grid information (one grid box, 256m thick, and 128x128 in width) and the computational architecture (CPU).

arch=CPU()

grid = RectilinearGrid(size=(1,1,1), x = (0, 128meters), y = (0,128meters), z = (0,-128meters))

# ## 3. Model Setup
# Next we setup the individual-based model by specifying the computational architecture, grid, and plankton community.

model = PlanktonModel(arch, grid; phyto = phyto_setup(nothing, [1024], 1),
                                  options = model_options(max_individuals = 2^10*16))

# And we setup diagnostics.

diags = PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                                   phytoplankton = (:num, :graz, :mort, :dvid, :PS, :BS, :Chl),
                                   iteration_interval = 1)

# Then we setup the duration of the model simulation, a run directory location, and the kind of output we want.

sim = PlanktonSimulation(model, ΔT = 300.0seconds,
                                iterations = 288,
                                diags = diags)

# Finally we setup the output writer.

writer = PlanktonOutputWriter(dir = joinpath(mktempdir(), "results"), save_diags = true)

# ## 4. Model Run

output_path = begin
    sim.output_writer = writer
    sim_run!(sim)
    writer.diags_file
end

# ## 5. Access Results
#
# Results have been stored in a `jld2` file. Let's open the file, look inside, and retrieve results.

series_names = jldopen(output_path, "r") do file
    keys(file["timeseries"])
end

# Extract a vector of iterations

iterations = jldopen(output_path, "r") do file
    sort!(parse.(Int, collect(keys(file["timeseries/t"]))))
end
nothing # hide

# Read results into arrays and close file.

function get_time_series!(file,name,array)
    f = jldopen(file, "r")
    for (i, iter) in enumerate(iterations)
        array[i]  = f["timeseries/$name/$iter"][1,1,1]
    end
    close(f)
end

begin
    (PAR, NH4, NO3, DOC) = (zeros(288),zeros(288),zeros(288),zeros(288))
    fil2 = output_path
    get_time_series!(fil2,"PAR",PAR)
    get_time_series!(fil2,"NH4",NH4)
    get_time_series!(fil2,"NO3",NO3)
    get_time_series!(fil2,"DOC",DOC)
end

begin
    (num,dvid,mort,PS,Chl) = (zeros(288),zeros(288),zeros(288),zeros(288),zeros(288))
    fil = output_path
    get_time_series!(fil,"phyto/sp1/num",num)
    get_time_series!(fil,"phyto/sp1/dvid",dvid)
    get_time_series!(fil,"phyto/sp1/mort",mort)
    get_time_series!(fil,"phyto/sp1/PS",PS)
    get_time_series!(fil,"phyto/sp1/Chl",Chl)
end

# ## 6. Visualize Results
# First we plot population size, division, mortality, and photosynthesis.

let
    p1 = plot(collect(1:300:300*288) ./ 3600, num, title = "population", legend=:none, fmt=:png, bottom_margin = 5mm)
    p2 = plot(collect(1:300:300*288) ./ 3600, dvid ./ num .* 60, title = "division rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
    p3 = plot(collect(1:300:300*288) ./ 3600, mort ./ num .* 60, title = "mortality rate (per hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
    p4 = plot(collect(1:300:300*288) ./ 3600, PS ./ num .* 12 .* 1e12 .* 3600, title = "photosynthesis rate (fg C/cell/hour)", legend=:none, fmt=:png, bottom_margin = 5mm)
    plot(p1, p2, p3, p4, layout = (4,1), size=(600,600), titlefont = (12))
end

# The second figure shows the light and nutrient conditions.

let
    p1 = plot(collect(1:300:300*288) ./ 3600, PAR, title = "PAR", legend=:none, fmt=:png, bottom_margin = 5mm)
    p2 = plot(collect(1:300:300*288) ./ 3600, NH4, title = "NH4", legend=:none, fmt=:png, bottom_margin = 5mm)
    p3 = plot(collect(1:300:300*288) ./ 3600, NO3, title = "NO3", legend=:none, fmt=:png, bottom_margin = 5mm)
    p4 = plot(collect(1:300:300*288) ./ 3600, DOC, title = "DOC", legend=:none, fmt=:png, bottom_margin = 5mm)
    p5 = plot(collect(1:300:300*288) ./ 3600, Chl, title = "Chl", legend=:none, fmt=:png, bottom_margin = 5mm)
    plot(p1, p2, p3, p4, p5, layout = (5,1), size=(600,600), titlefont = (12))
end
