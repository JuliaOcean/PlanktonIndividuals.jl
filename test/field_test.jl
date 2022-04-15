using PlanktonIndividuals.Fields
using PlanktonIndividuals.Biogeochemistry


function test_fields()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))

    nut = nutrients_init(CPU(), grid)

    @test Tuple(collect(keys(nut))) == nut_names
    @test nut.DIC.data == zeros(8,10,6)
    @test interior(nut.DIC.data, grid) == zeros(4,6,2) 
    
    nuts = generate_nutrients(CPU(), grid, default_nut_init())
    @test maximum(nuts.DIC.data) < 23.0
    @test minimum(nuts.DIC.data) > 17.0

    zero_fields!(nuts)
    @test nuts.DIC.data == zeros(8,10,6)

    return nothing
end
function test_fill_halos()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8))
    nuts = generate_nutrients(CPU(), grid, default_nut_init())
    Nx,Ny,Nz = grid.Nx, grid.Ny, grid.Nz
    Hx,Hy,Hz = grid.Hx, grid.Hy, grid.Hz

    PlanktonIndividuals.Fields.fill_halo_east!(nuts.DIC.data, Hx, Nx, Periodic())
    @test nuts.DIC.data[Nx+Hx+1:Nx+2*Hx, :, :] == nuts.DIC.data[1+Hx:2*Hx, :, :]

    PlanktonIndividuals.Fields.fill_halo_west!(nuts.DIC.data, Hx, Nx, Periodic())
    @test nuts.DIC.data[1:Hx, :, :] == nuts.DIC.data[Nx+1:Nx+Hx, :, :]

    PlanktonIndividuals.Fields.fill_halo_north!(nuts.DIC.data, Hy, Ny, Periodic())
    @test nuts.DIC.data[:, Ny+Hy+1:Ny+2*Hy, :] == nuts.DIC.data[:, 1+Hy:2*Hy, :]

    PlanktonIndividuals.Fields.fill_halo_south!(nuts.DIC.data, Hy, Ny, Periodic())
    @test nuts.DIC.data[:, 1:Hy, :] == nuts.DIC.data[:, Ny+1:Ny+Hy, :]

    PlanktonIndividuals.Fields.fill_halo_bottom!(nuts.DIC.data, Hz, Nz, Periodic())
    @test nuts.DIC.data[:, :, Nz+Hz+1:Nz+2*Hz] == nuts.DIC.data[:, :, 1+Hz:2*Hz]

    PlanktonIndividuals.Fields.fill_halo_top!(nuts.DIC.data, Hz, Nz, Periodic())
    @test nuts.DIC.data[:, :, 1:Hz] == nuts.DIC.data[:, :, Nz+1:Nz+Hz]

    PlanktonIndividuals.Fields.fill_halo_east!(nuts.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test nuts.DIC.data[Nx+Hx+i, :, :] == nuts.DIC.data[Nx+Hx, :, :]
    end

    PlanktonIndividuals.Fields.fill_halo_west!(nuts.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test nuts.DIC.data[i, :, :] == nuts.DIC.data[Hx+1, :, :]
    end

    PlanktonIndividuals.Fields.fill_halo_north!(nuts.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test nuts.DIC.data[:, Ny+Hy+i, :] == nuts.DIC.data[:, Ny+Hy, :]
    end

    PlanktonIndividuals.Fields.fill_halo_south!(nuts.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test nuts.DIC.data[:, i, :] == nuts.DIC.data[:, Hy+1, :]
    end

    PlanktonIndividuals.Fields.fill_halo_bottom!(nuts.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test nuts.DIC.data[:, :, Nz+Hz+i] == nuts.DIC.data[:, :, Nz+Hz]
    end

    PlanktonIndividuals.Fields.fill_halo_top!(nuts.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test nuts.DIC.data[:, :, i] == nuts.DIC.data[:, :, Hz+1]
    end

    PlanktonIndividuals.Fields.fill_halo_east_vel!(nuts.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx-1
        @test nuts.DIC.data[Nx+Hx+1+i, :, :] == nuts.DIC.data[Nx+Hx+1, :, :]
    end

    PlanktonIndividuals.Fields.fill_halo_north_vel!(nuts.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy-1
        @test nuts.DIC.data[:, Ny+Hx+1+i, :] == nuts.DIC.data[:, Ny+Hy+1, :]
    end

    PlanktonIndividuals.Fields.fill_halo_bottom_vel!(nuts.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz-1
        @test nuts.DIC.data[:, :, Nz+Hx+1+i] == nuts.DIC.data[:, :, Nz+Hz+1]
    end

    PlanktonIndividuals.Fields.fill_halo_east_Gc!(nuts.DIC.data, Hx, Nx, Bounded())
    for i in 1:Hx
        @test nuts.DIC.data[Nx+Hx+i, :, :] == zeros(10,6)
    end

    PlanktonIndividuals.Fields.fill_halo_north_Gc!(nuts.DIC.data, Hy, Ny, Bounded())
    for i in 1:Hy
        @test nuts.DIC.data[:, Ny+Hx+i, :] == zeros(8,6)
    end

    PlanktonIndividuals.Fields.fill_halo_bottom_Gc!(nuts.DIC.data, Hz, Nz, Bounded())
    for i in 1:Hz
        @test nuts.DIC.data[:, :, Nz+Hx+i] == zeros(8,10)
    end

    return nothing
end

function test_boundary_conditions()
    grid = RectilinearGrid(size = (4, 4, 4), x = (0,32), y = (0,32), z = (0,-32))
    model = PlanktonModel(CPU(), grid; mode = CarbonMode()) 
    set_bc!(model, :DIC, :west, 0.1)
    @test model.nutrients.DIC.bc.west == 0.1
    set_bc!(model, :DIC, :west, ones(4,4))
    @test model.nutrients.DIC.bc.west == ones(4,4)
    set_bc!(model, :DIC, :west, ones(4,4,10))
    @test model.nutrients.DIC.bc.west == ones(4,4,10)

    Gcs = nutrients_init(CPU(), grid)
    apply_bcs!(Gcs, model.nutrients, model.grid, 10, 1, CPU())
    @test Gcs.DIC.data[3,3:6,3:6] == ones(4,4) ./ 8.0

end

@testset "Fields" begin
    @testset "Fields" begin
        test_fields()
    end
    @testset "Fill halos" begin
        test_fill_halos()
    end
    @testset "Boundary conditions" begin
        test_boundary_conditions()
    end
end