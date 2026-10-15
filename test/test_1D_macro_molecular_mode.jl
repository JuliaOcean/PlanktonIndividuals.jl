@testset "PlanktonIndividuals 1D MacroMolecularMode tests" begin
    Random.seed!(1234)
    grid = RectilinearGrid(size = (1, 1, 16), x = (0,32), y = (0,32), z = (0,-32),
                                  topology = (Bounded, Bounded, Bounded), halo = (2, 2, 2))

    model = PlanktonModel(CPU(), grid; mode = MacroMolecularMode(),
                         phyto = phyto_setup(Dict("κvP" => [1e-4]), [1024], 1),
                         bgc_params = Dict("κv" => 1e-4))

    uvel = zeros(2,1,16,11)
    vvel = zeros(1,2,16,11)
    wvel = zeros(1,1,17,11)

    # Bounded column: impermeable faces; particle/tracer diffusion and physiology remain active.

    test_example_state(model; initial_counts=[1024])

    sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10, vels=(u=uvel, v=vvel, w=wvel))

    test_example_evolution!(sim)
end
