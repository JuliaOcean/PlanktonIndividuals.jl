using PlanktonIndividuals, PlanktonKernels
using PlanktonKernels.Fields: BoundaryCondition, BoundaryConditions

@testset "Abiotic boundary fluxes" begin
    ab = PlanktonIndividuals.Individuals.Abiotic
    for FT in (Float32, Float64)
        g = RectilinearGrid(size=(3,4,2), x=(0,6), y=(0,12), z=(0,-8), FT=FT)
        dims = (g.Nx+2g.Hx,g.Ny+2g.Hy,g.Nz+2g.Hz)
        p = (sz_min=FT(2), Nsuper=FT(3))
        rnd = fill(FT(0.5),dims)
        for (face,axis,upper,shape,area) in
            ((:west,1,false,(4,2),12),(:east,1,true,(4,2),12),
             (:south,2,false,(3,2),8),(:north,2,true,(3,2),8),
             (:top,3,false,(3,4),6),(:bottom,3,true,(3,4),6))
            @testset "$FT $face" begin
                spatial = reshape(FT.(1:prod(shape)),shape) ./ FT(10)
                temporal = cat(spatial, reverse(spatial); dims=3)
                for c in (zero(FT),FT(-1),FT(0.5),spatial,temporal),
                    v in (zero(FT),FT(-1),FT(0.5),spatial,temporal), t in 1:2
                    # Compare each cell against a scalar probability threshold.
                    out = fill(FT(-7),dims)
                    expected = copy(out)
                    for b in 1:shape[2], a in 1:shape[1]
                        cv = c isa Number ? c : ndims(c)==2 ? c[a,b] : c[a,b,t]
                        vv = v isa Number ? v : ndims(v)==2 ? v[a,b] : v[a,b,t]
                        edge = upper ? (g.Nx,g.Ny,g.Nz)[axis] : 1
                        idx = axis==1 ? (edge,a,b) : axis==2 ? (a,edge,b) : (a,b,edge)
                        ii,jj,kk = idx .+ (g.Hx,g.Hy,g.Hz)
                        expected[ii,jj,kk] = FT(3) < cv*vv*FT(0.5)*FT(area)
                    end
                    bcs = BoundaryConditions(nothing,nothing,nothing,nothing,nothing,nothing)
                    setproperty!(bcs,face,BoundaryCondition(c,v))
                    ab.calc_particle_bcs!(out,bcs,rnd,p,FT(0.5),t,g,CPU())
                    @test out == expected
                end
            end
        end
        out = fill(FT(7),dims)
        bcs = BoundaryConditions(nothing,nothing,nothing,nothing,nothing,nothing)
        ab.calc_particle_bcs!(out,bcs,rnd,p,FT(0.5),1,g,CPU())
        @test all(==(FT(7)),out)
    end
end

@testset "Abiotic boundary setter" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            bcs = BoundaryConditions(nothing,nothing,nothing,nothing,nothing,nothing)
            particle = PlanktonIndividuals.AbioticParticle(zeros(FT,1), (;), bcs)
            @test PlanktonIndividuals.Individuals.Abiotic.set_bc_particle!(particle, :top, FT(0.5), FT, CPU()) === nothing
            @test particle.bc.top isa BoundaryCondition
        end
    end
end
