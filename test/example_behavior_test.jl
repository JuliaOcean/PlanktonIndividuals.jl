# Controlled experiments complement the larger examples: use identical initial
# cells for each comparison and suppress demographic events to isolate physiology.
function physiology_example(mode; FT=Float32, counts=[16], size=(1, 1, 1))
    # Give uptake a resolvable signal against dissolved background in Float32.
    params = Dict("Nsuper" => fill(1e10, length(counts)), "grz_P" => zeros(length(counts)),
                  "mort_P" => zeros(length(counts)),
                  "dvid_P" => zeros(length(counts)))
    if mode isa ProteinMode
        # Isolate uptake from maintenance-energy starvation. The larger default
        # ProteinMode example above still exercises its normal maintenance cost.
        params["e_min"] = zeros(length(counts))
    end
    grid = RectilinearGrid(; size, x=(0, 2), y=(0, 2), z=(0, -2))
    return PlanktonModel(CPU(), grid; mode, FT,
                         phyto=phyto_setup(params, counts, length(counts)),
                         options=model_options(max_individuals=64),
                         bgc_params=Dict("w_sink_org"=>0.0, "w_sink_inorg"=>0.0))
end

function forced_example(model, light; iterations=1, vels=(;))
    g = model.grid
    return PlanktonSimulation(model; ΔT=10.0, iterations, vels,
                              PARF=fill(light, g.Nx, g.Ny, 24),
                              temp=fill(20.0, g.Nx, g.Ny, g.Nz, 24))
end

@testset "Controlled example physiology" begin
    for mode in (CarbonMode(), QuotaMode(), MacroMolecularMode(), IronEnergyMode(), ProteinMode())
        @testset "$(typeof(mode))" begin
            Random.seed!(2026)
            initial = physiology_example(mode)
            @test typeof(initial.mode) === typeof(mode)
            pools = mode isa CarbonMode ? (:Bm, :Bd) :
                    mode isa QuotaMode ? (:Cq, :Nq, :Pq) :
                    mode isa MacroMolecularMode ? (:PRO, :DNA, :RNA, :NST, :PST) :
                    mode isa IronEnergyMode ? (:ATP, :NADPH, :qFePS, :qFeNR, :qFeNF) :
                    (:PRO_RB, :PRO_PS, :PRO_MC, :Chl_lag)
            @test all(key -> hasproperty(initial.individuals.phytos.sp1.data, key), pools)
            test_example_state(initial; initial_counts=[16])
            for FT in (Float32, Float64)
                # Exercise precision dispatch separately; compare directions, not trajectories.
                Random.seed!(2026)
                dark = physiology_example(mode; FT)
                lit = deepcopy(dark)
                dark_sim, lit_sim = forced_example(dark, 0.0), forced_example(lit, 200.0)
                initial_inventory = example_inventory(dark)
                initial_P = mode isa CarbonMode ? 0.0 : particle_phosphorus(dark.individuals.phytos.sp1, mode)
                update!(dark_sim)
                update!(lit_sim)
                for model in (dark, lit)
                    d = model.individuals.phytos.sp1.data
                    @test count(d.ac) == 16
                    @test all(==(1), d.gen[d.ac])
                    @test all(x -> isapprox(x, 10/3600; rtol=1e-6), d.age[d.ac])
                    test_example_state(model)
                    @test example_inventory(model) ≈ initial_inventory rtol=5e-6
                end
                dd, ld = dark.individuals.phytos.sp1.data, lit.individuals.phytos.sp1.data
                @test all(iszero, dd.PS[dd.ac])
                @test all(>(0), ld.PS[ld.ac])
                @test all(iszero, ld.PS[.!ld.ac])
                @test all(iszero, dark.timestepper.par)
                @test maximum(lit.timestepper.par) > 0
                if mode isa CarbonMode
                    # Fixed Chl:C is a design feature only of CarbonMode.
                    for model in (dark, lit)
                        d, p = model.individuals.phytos.sp1.data, model.individuals.phytos.sp1.p
                        @test d.Chl ≈ d.Bm .* p.Chl2C
                        @test d.Bm ≈ d.Sz .* (p.Cquota * p.Nsuper)
                        @test all(d.Bd .<= d.Bm)
                    end
                    @test sum(ld.Bm) > sum(dd.Bm)
                    # More fixation removes more DIC with all other processes matched.
                    @test tot_mass(lit.tracers.DIC.data, lit.grid) < tot_mass(dark.tracers.DIC.data, dark.grid)
                end
                if !(mode isa CarbonMode)
                    @test sum(ld.VPO4[ld.ac]) > 0
                    final_P = particle_phosphorus(lit.individuals.phytos.sp1, mode)
                    @test final_P > initial_P
                    if FT === Float64
                        # Check the transfer itself, so a large dissolved background
                        # cannot conceal missing particle uptake or field coupling.
                        uptake = 10 * sum(ld.VPO4[ld.ac])
                        @test final_P - initial_P ≈ uptake rtol=1e-4
                        dissolved_P = sum(tot_mass(lit.tracers[key].data, lit.grid)
                                          for key in (:PO4, :DOP, :POP))
                        @test dissolved_P - (initial_inventory - initial_P) ≈ -uptake rtol=1e-4
                    end
                end
                if mode isa IronEnergyMode || mode isa ProteinMode
                    # Defaults disable diazotrophy: its process rate must remain zero.
                    p = lit.individuals.phytos.sp1.p
                    if p.is_croc + p.is_tric == 0
                        @test all(iszero, ld.NF[ld.ac])
                    end
                end
            end
        end
    end
end

@testset "Equivalent example configurations" begin
    for mode in (CarbonMode(), QuotaMode(), MacroMolecularMode(), IronEnergyMode(), ProteinMode())
        @testset "$(typeof(mode))" begin
            Random.seed!(42)
            model = physiology_example(mode; counts=[8, 12], size=(2, 2, 2))
            test_example_state(model; initial_counts=[8, 12])
            batch, split, zero_flow = deepcopy(model), deepcopy(model), deepcopy(model)
            g = model.grid
            # Periodic x/y and bounded z require Nz+1 vertical faces.
            vels = (u=zeros(g.Nx, g.Ny, g.Nz, 3),
                    v=zeros(g.Nx, g.Ny, g.Nz, 3),
                    w=zeros(g.Nx, g.Ny, g.Nz+1, 3))
            update!(forced_example(batch, 200.0; iterations=2))
            split_sim = forced_example(split, 200.0)
            update!(split_sim)
            update!(split_sim)
            update!(forced_example(zero_flow, 200.0; iterations=2, vels))
            for other in (split, zero_flow)
                @test other.t == batch.t == 20
                @test other.iteration == batch.iteration == 2
                for key in keys(batch.tracers)
                    @test other.tracers[key].data ≈ batch.tracers[key].data
                end
                for key in keys(batch.individuals.phytos)
                    d, ref = other.individuals.phytos[key].data, batch.individuals.phytos[key].data
                    @test d.ac == ref.ac
                    @test d.Chl ≈ ref.Chl
                    @test d.PS[d.ac] ≈ ref.PS[ref.ac]
                    @test d.x ≈ ref.x
                    @test d.y ≈ ref.y
                    @test d.z ≈ ref.z
                end
            end
            # Verify exchange for every species in a small-volume domain, where a
            # missing or repeated amount-to-concentration conversion is visible.
            @test example_inventory(batch) ≈ example_inventory(model) rtol=5e-6
        end
    end
end

@testset "Quota nutrient limitation and column light" begin
    model = physiology_example(QuotaMode(); size=(1, 1, 4))
    depleted = deepcopy(model)
    depleted.tracers.PO4.data .= 0
    update!(forced_example(model, 200.0))
    update!(forced_example(depleted, 200.0))
    rich, poor = model.individuals.phytos.sp1.data, depleted.individuals.phytos.sp1.data
    @test all(>(0), rich.VPO4[rich.ac])
    @test all(iszero, poor.VPO4[poor.ac])
    # Cell-averaged PAR is positive, less than incident PAR, and decreases with depth.
    g = model.grid
    column = model.timestepper.par[g.Hx+1, g.Hy+1, g.Hz+1:g.Hz+g.Nz]
    @test all(x -> 0 < x < 200, column)
    @test all(<(0), diff(column))
end

@testset "Example boundary flux mass balance" begin
    grid = RectilinearGrid(size=(2, 3, 4), x=(0, 2), y=(0, 6), z=(0, -8),
                           topology=(Bounded, Bounded, Bounded))
    model = PlanktonModel(CPU(), grid; mode=CarbonMode(), FT=Float64,
                         phyto=phyto_setup(Dict("grz_P"=>[0.0], "mort_P"=>[0.0], "dvid_P"=>[0.0]), [8], 1),
                         options=model_options(max_individuals=32),
                         bgc_params=Dict("kDON"=>0.0, "w_sink_org"=>0.0, "w_sink_inorg"=>0.0))
    # Retain scalar, spatial, and time-dependent boundary coverage from the 3D
    # example, now with an analytic budget and no DON reaction in CarbonMode.
    model.tracers.DON.bc.west = 1e-3
    model.tracers.DON.bc.east = fill(2e-4, 3, 4)
    model.tracers.DON.bc.south = cat(fill(3e-4, 2, 4), fill(6e-4, 2, 4); dims=3)
    model.tracers.DON.bc.north = cat(fill(1e-4, 2, 4), fill(2e-4, 2, 4); dims=3)
    initial_DON = tot_mass(model.tracers.DON.data, model.grid)
    initial_C = example_inventory(model)
    update!(forced_example(model, 0.0; iterations=2))
    # Signed inflow × face area × duration. Both time slices must be applied once.
    expected = (1e-3 - 2e-4) * (6 * 8) * 20 +
               ((3e-4 - 1e-4) + (6e-4 - 2e-4)) * (2 * 8) * 10
    @test tot_mass(model.tracers.DON.data, model.grid) - initial_DON ≈ expected rtol=1e-6
    @test example_inventory(model) ≈ initial_C rtol=5e-6
end
