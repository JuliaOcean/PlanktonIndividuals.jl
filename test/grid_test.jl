using PlanktonIndividuals.Grids

using MeshArrays

function test_rectilinear_grid()
    grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8), halo = (2,2,2))

    @test grid.Nx == 4
    @test grid.Ny == 6
    @test grid.Nz == 2

    @test grid.Δx == 3.0
    @test grid.Δy == 2.0
    @test grid.dzC == [4.0,4.0,4.0,4.0,4.0,4.0]
    @test grid.dzF == [4.0,4.0,4.0,4.0,4.0,4.0]

    @test length(grid.xC) == 4+2*2
    @test length(grid.yC) == 6+2*2
    @test length(grid.zC) == 2+2*2
    @test length(grid.xF) == 4+2*2
    @test length(grid.yF) == 6+2*2
    @test length(grid.zF) == 2+2*2

    return nothing
end

function test_vertically_stretched_rectilinear_grid()
    zfs = [0.0, -1.0, -3.0, -5.0, -10.0]
    grid = RectilinearGrid(size = (4,6,4), x = (0,12), y = (0,12), z = zfs, halo = (2,2,2))

    @test grid.Nx == 4
    @test grid.Ny == 6
    @test grid.Nz == 4

    @test grid.Δx == 3.0
    @test grid.Δy == 2.0
    @test grid.dzC == [1.0,1.0,1.5,2.0,3.5,5.0,5.0,5.0]
    @test grid.dzF == [1.0,1.0,1.0,2.0,2.0,5.0,5.0,5.0]

    @test grid.zF == [2.0,1.0,-0.0,-1.0,-3.0,-5.0,-10.0,-15.0]
    @test grid.zC == [1.5,0.5,-0.5,-2.0,-4.0,-7.5,-12.5,-17.5]

    @test length(grid.xC) == 4+2*2
    @test length(grid.yC) == 6+2*2
    @test length(grid.zC) == 4+2*2
    @test length(grid.xF) == 4+2*2
    @test length(grid.yF) == 6+2*2
    @test length(grid.zF) == 4+2*2

    return nothing
end

function test_rectilinear_areas_volumes()
    zfs = [0.0, -1.0, -3.0, -5.0, -10.0]
    grid = RectilinearGrid(size = (4,6,4), x = (0,12), y = (0,12), z = zfs, halo = (2,2,2))
    @test ΔxC(1,1,1,grid) == grid.Δx
    @test ΔyC(1,1,1,grid) == grid.Δy
    @test ΔzC(1,1,1,grid) == grid.dzC[1]
    @test ΔxF(1,1,1,grid) == grid.Δx
    @test ΔyF(1,1,1,grid) == grid.Δy
    @test ΔzF(1,1,1,grid) == grid.dzF[1]
    @test Ax(1,1,1,grid) == grid.Δy*grid.dzF[1]
    @test Ay(1,1,1,grid) == grid.Δx*grid.dzF[1]
    @test Az(1,1,1,grid) == grid.Δx*grid.Δy
    @test volume(1,1,1,grid) == grid.Δx*grid.Δy*grid.dzF[1]

    return nothing
end

function test_lat_lon_grid()
    grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
    @test grid.Nx == 360
    @test grid.Ny == 160
    @test grid.Nz == 10

    @test length(grid.xC) == 360+2*2
    @test length(grid.yC) == 160+2*2
    @test length(grid.zC) == 10+2*2
    @test length(grid.xF) == 360+2*2
    @test length(grid.yF) == 160+2*2
    @test length(grid.zF) == 10+2*2
end

function test_load_lat_lon_grid()
    MeshArrays.GRID_LL360_download();
	γ = MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360);
	Γ = MeshArrays.GridLoad(γ;option="full");
	
	grid_info=(XC=Γ.XC[1], YC=Γ.YC[1], RC=Γ.RC,
        	XW=Γ.XW[1], YS=Γ.YS[1], RF=Γ.RF,
            DXC=Γ.DXC[1], DYC=Γ.DYC[1], DRC=Γ.DRC,
            DXG=Γ.DXG[1], DYG=Γ.DYG[1], DRF=Γ.DRF,
            RAW=Γ.RAW[1], RAS=Γ.RAS[1], RAC=Γ.RAC[1],
            hFacC=ones(360,160,50), hFacW=ones(360,160,50), hFacS=ones(360,160,50))

    mask = ones(360,160,50)

    grid = LoadLatLonGrid(;grid_info = grid_info, size = (360,160,50), lat = (-80,80), lon= (-180,180), landmask = mask)

    @test grid.Nx == 360
    @test grid.Ny == 160
    @test grid.Nz == 50

    @test length(grid.xC) == 360+2*2
    @test length(grid.yC) == 160+2*2
    @test length(grid.zC) == 50+2*2
    @test length(grid.xF) == 360+2*2
    @test length(grid.yF) == 160+2*2
    @test length(grid.zF) == 50+2*2

    @test minimum(grid.dzF) == 10.0
    @test maximum(grid.dzF) == 456.5

    return nothing
end

function test_lat_lon_areas_volumes()
    grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
    @test ΔxC(1,1,1,grid) == grid.dxC[1,1]
    @test ΔyC(1,1,1,grid) == grid.dyC[1,1]
    @test ΔzC(1,1,1,grid) == grid.dzC[1]
    @test ΔxF(1,1,1,grid) == grid.dxF[1,1]
    @test ΔyF(1,1,1,grid) == grid.dyF[1,1]
    @test ΔzF(1,1,1,grid) == grid.dzF[1]
    @test Ax(1,1,1,grid) == grid.Ax[1,1,1]
    @test Ay(1,1,1,grid) == grid.Ay[1,1,1]
    @test Az(1,1,1,grid) == grid.Az[1,1]
    @test volume(1,1,1,grid) == grid.Vol[1,1,1]

    return nothing
end

@testset "Grids" begin
    @testset "Rectilinear Grid" begin
        test_rectilinear_grid()
        test_vertically_stretched_rectilinear_grid()
        test_rectilinear_areas_volumes()
        grid = RectilinearGrid(size = (4,6,2), x = (0,12), y = (0,12), z = (0,-8), halo = (2,2,2))
        @test try
            show(grid); println()
            true
        catch err
            println("error in show(::RectilinearGrid)")
            println(sprint(showerror, err))
            false
        end
        @test grid isa RectilinearGrid
    end
    @testset "Latitude Longitude Grid" begin
        test_lat_lon_grid()
        test_load_lat_lon_grid()
        test_lat_lon_areas_volumes()
        grid = LatLonGrid(size = (360,160,10), lat = (-80,80), lon = (-180,180), z = (0,-20))
        @test try
            show(grid); println()
            true
        catch err
            println("error in show(::LatLonGrid)")
            println(sprint(showerror, err))
            false
        end
        @test grid isa LatLonGrid
    end
end