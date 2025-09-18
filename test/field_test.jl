using PlanktonIndividuals.Biogeochemistry

function test_fields()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))

    tracers = tracers_init(CPU(), grid)

    @test Tuple(collect(keys(tracers))) == tracer_names
    @test tracers.DIC.data == zeros(8,10,6)
    @test interior(tracers.DIC.data, grid) == zeros(4,6,2) 
    
    tracers = generate_tracers(CPU(), grid, default_tracer_init(), Float32)
    @test maximum(tracers.DIC.data) < 23.0
    @test minimum(tracers.DIC.data) > 17.0

    zero_fields!(tracers)
    @test tracers.DIC.data == zeros(8,10,6)

    return nothing
end
function test_fill_halos()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))
    tracers = generate_tracers(CPU(), grid, default_tracer_init(), Float32)
    Nx,Ny,Nz = grid.Nx, grid.Ny, grid.Nz
    Hx,Hy,Hz = grid.Hx, grid.Hy, grid.Hz

    PlanktonIndividuals.Biogeochemistry.fill_halo_east!(tracers.DIC.data, Hx, Nx, Periodic())
    @test tracers.DIC.data[Nx+Hx+1:Nx+2*Hx, :, :] == tracers.DIC.data[1+Hx:2*Hx, :, :]

    PlanktonIndividuals.Biogeochemistry.fill_halo_west!(tracers.DIC.data, Hx, Nx, Periodic())
    @test tracers.DIC.data[1:Hx, :, :] == tracers.DIC.data[Nx+1:Nx+Hx, :, :]

    PlanktonIndividuals.Biogeochemistry.fill_halo_north!(tracers.DIC.data, Hy, Ny, Periodic())
    @test tracers.DIC.data[:, Ny+Hy+1:Ny+2*Hy, :] == tracers.DIC.data[:, 1+Hy:2*Hy, :]

    PlanktonIndividuals.Biogeochemistry.fill_halo_south!(tracers.DIC.data, Hy, Ny, Periodic())
    @test tracers.DIC.data[:, 1:Hy, :] == tracers.DIC.data[:, Ny+1:Ny+Hy, :]

    PlanktonIndividuals.Biogeochemistry.fill_halo_bottom!(tracers.DIC.data, Hz, Nz, Periodic())
    @test tracers.DIC.data[:, :, Nz+Hz+1:Nz+2*Hz] == tracers.DIC.data[:, :, 1+Hz:2*Hz]

    PlanktonIndividuals.Biogeochemistry.fill_halo_top!(tracers.DIC.data, Hz, Nz, Periodic())
    @test tracers.DIC.data[:, :, 1:Hz] == tracers.DIC.data[:, :, Nz+1:Nz+Hz]

    PlanktonIndividuals.Biogeochemistry.fill_halo_east!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[Nx+Hx+i, :, :] == tracers.DIC.data[Nx+Hx, :, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_west!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[i, :, :] == tracers.DIC.data[Hx+1, :, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_north!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, Ny+Hy+i, :] == tracers.DIC.data[:, Ny+Hy, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_south!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, i, :] == tracers.DIC.data[:, Hy+1, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_bottom!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, Nz+Hz+i] == tracers.DIC.data[:, :, Nz+Hz]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_top!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, i] == tracers.DIC.data[:, :, Hz+1]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_east_vel!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx-1
        @test tracers.DIC.data[Nx+Hx+1+i, :, :] == tracers.DIC.data[Nx+Hx+1, :, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_north_vel!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy-1
        @test tracers.DIC.data[:, Ny+Hx+1+i, :] == tracers.DIC.data[:, Ny+Hy+1, :]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_bottom_vel!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz-1
        @test tracers.DIC.data[:, :, Nz+Hx+1+i] == tracers.DIC.data[:, :, Nz+Hz+1]
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_east_Gc!(tracers.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test tracers.DIC.data[Nx+Hx+i, :, :] == zeros(10,6)
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_north_Gc!(tracers.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test tracers.DIC.data[:, Ny+Hx+i, :] == zeros(8,6)
    end

    PlanktonIndividuals.Biogeochemistry.fill_halo_bottom_Gc!(tracers.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test tracers.DIC.data[:, :, Nz+Hx+i] == zeros(8,10)
    end

    return nothing
end

function test_boundary_conditions()
    grid = RectilinearGrid(size = (4, 4, 4), x = (0,32), y = (0,32), z = (0,-32))
    model = PlanktonModel(CPU(), grid; mode = CarbonMode()) 
    FT = model.FT
    set_bc!(model; tracer = :DIC, pos = :west, bc_value = 0.1)
    @test model.tracers.DIC.bc.west == FT(0.1)
    set_bc!(model; tracer = :DIC, pos = :west, bc_value = ones(4,4))
    @test model.tracers.DIC.bc.west == ones(FT, 4,4)
    set_bc!(model; tracer = :DIC, pos = :west, bc_value = ones(4,4,10))
    @test model.tracers.DIC.bc.west == ones(FT, 4,4,10)

    Gcs = tracers_init(CPU(), grid)
    apply_bcs!(Gcs, model.tracers, model.grid, 10, 1, CPU())
    @test Gcs.DIC.data[3,3:6,3:6] == ones(FT, 4,4) ./ FT(8.0)

end

@testset "Biogeochemistry" begin
    test_fields()
    test_fill_halos()
    test_boundary_conditions() 
end