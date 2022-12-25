using PlanktonIndividuals.Output
using JLD2

function test_output()
    grid = RectilinearGrid(size = (16, 16, 16), x = (0,32), y = (0,32), z = (0,-32))
    model = PlanktonModel(CPU(), grid; mode = QuotaMode()) 
    diags = PlanktonDiagnostics(model; tracer=(:PAR,),
                                   plankton = (:num, :graz, :mort, :dvid, :PS),
                                   iteration_interval = 1)

    sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 4, diags = diags) 

    sim.output_writer = PlanktonOutputWriter(dir = "./result",
                                             write_log = true,
                                             save_diags = true,
                                             save_plankton = true,
                                             plankton_iteration_interval = 1,
                                             max_filesize = 256KiB)
    
    @test sim.output_writer isa PlanktonOutputWriter

    update!(sim)
    
    @test isfile("result/dynamic_species001.txt")
    @test isfile("result/diags_part1.jld2")
    @test isfile("result/plankton_part1.jld2")

    ds = jldopen("result/diags_part1.jld2")
    time_length = length(keys(ds["timeseries/t"]))
    @test time_length == 2
    close(ds)

    ds1 = jldopen("result/plankton_part1.jld2")
    time_length1 = length(keys(ds1["timeseries/t"]))
    @test time_length1 == 1
    close(ds1)

    rm("result", recursive=true)
    return nothing
end

@testset "Output" begin
    test_output()
end
