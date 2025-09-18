using PlanktonIndividuals.Output
using JLD2

function test_output()
    grid = RectilinearGrid(size = (16, 16, 16), x = (0,32), y = (0,32), z = (0,-32))
    model = PlanktonModel(CPU(), grid; mode = QuotaMode(),
                                       abiotic = (params = nothing, N = [2^10], Nsa = 1)) 
    diags = PlanktonDiagnostics(model; tracer=(:PAR,),
                                       phytoplankton = (:num, :graz, :mort, :dvid, :PS),
                                       abiotic_particle = (:num, :CHO),
                                       iteration_interval = 1)

    sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 4, diags = diags) 

    sim.output_writer = PlanktonOutputWriter(dir = "./result",
                                             write_log = true,
                                             save_diags = true,
                                             save_phytoplankton = true,
                                             save_abiotic_particle = true,
                                             phytoplankton_iteration_interval = 1,
                                             abiotic_particle_iteration_interval = 1,
                                             max_filesize = 256KiB)
    
    @test sim.output_writer isa PlanktonOutputWriter

    update!(sim)
    
    @test isfile("result/dynamic_species001.txt")
    @test isfile("result/diags_part1.jld2")
    @test isfile("result/phytoplankton_part1.jld2")
    @test isfile("result/abiotic_particle_part1.jld2")

    ds = jldopen("result/diags_part1.jld2")
    time_length = length(keys(ds["timeseries/t"]))
    @test time_length > 1
    close(ds)

    ds1 = jldopen("result/phytoplankton_part1.jld2")
    time_length1 = length(keys(ds1["timeseries/t"]))
    @test time_length1 > 1
    close(ds1)

    ds2 = jldopen("result/abiotic_particle_part1.jld2")
    time_length2 = length(keys(ds2["timeseries/t"]))
    @test time_length2 > 1
    close(ds2)

    rm("result", recursive=true)
    return nothing
end

@testset "Output" begin
    test_output()
end
