@testset "PlanktonIndividuals 2D QuotaMode tests" begin
    Random.seed!(1234)
    grid = RectilinearGrid(size = (16, 16, 1), x = (0,32), y = (0,32), z = (0,-32))

    model = PlanktonModel(CPU(), grid; abiotic = abiotic_setup(nothing, [2^10, 2^10], 2, Palat([(:sp1, :sa1)], [(:sp1, :sa1)])))

    uvel = zeros(16,16,1,11)
    vvel = zeros(16,16,1,11)
    wvel = zeros(16,16,2,11)

    # Uniform horizontal flow is divergence-free on the periodic domain.
    uvel .= 1e-4
    vvel .= -1e-4

    test_example_state(model; initial_counts=[1024])

    sim = PlanktonSimulation(model, ΔT = 60.0, iterations = 10, vels=(u=uvel, v=vvel, w=wvel))

    test_example_evolution!(sim)
end
