@testset "PlanktonIndividuals 0D ProteinMode tests" begin
    Random.seed!(1234)
    grid = RectilinearGrid(size = (1, 1, 1), x = (0,32), y = (0,32), z = (0,-32))

    model = PlanktonModel(CPU(), grid;
                          mode = ProteinMode(),
                          phyto = phyto_setup(nothing, [1024], 1),
                          options = let opt = model_options(); opt.max_individuals = 1024*10; opt end)

    test_example_state(model; initial_counts=[1024])

    sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10)

    test_example_evolution!(sim)
end
