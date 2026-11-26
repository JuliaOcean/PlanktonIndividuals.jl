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
                sim_run!(dark_sim)
                sim_run!(lit_sim)
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
            sim_run!(forced_example(batch, 200.0; iterations=2))
            split_sim = forced_example(split, 200.0)
            sim_run!(split_sim)
            sim_run!(split_sim)
            sim_run!(forced_example(zero_flow, 200.0; iterations=2, vels))
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
    sim_run!(forced_example(model, 200.0))
    sim_run!(forced_example(depleted, 200.0))
    rich, poor = model.individuals.phytos.sp1.data, depleted.individuals.phytos.sp1.data
    @test all(>(0), rich.VPO4[rich.ac])
    @test all(iszero, poor.VPO4[poor.ac])
    # Cell-averaged PAR is positive, less than incident PAR, and decreases with depth.
    g = model.grid
    column = model.timestepper.par[g.Hx+1, g.Hy+1, g.Hz+1:g.Hz+g.Nz]
    @test all(x -> 0 < x < 200, column)
    @test all(<(0), diff(column))
end

function carbon_mode_nutrient_inventory(model)
    phyto = model.individuals.phytos.sp1
    biomass = sum(phyto.data.Bm[phyto.data.ac])
    nitrogen = sum(tot_mass(model.tracers[key].data, model.grid) for key in (:NH4, :NO3, :DON, :PON))
    phosphorus = sum(tot_mass(model.tracers[key].data, model.grid) for key in (:PO4, :DOP, :POP))
    return [nitrogen + biomass * phyto.p.R_NC, phosphorus + biomass * phyto.p.R_PC]
end

@testset "CarbonMode nutrient limitation" begin
    Random.seed!(2026)
    initial = physiology_example(CarbonMode(); FT=Float64)
    phyto = initial.individuals.phytos.sp1
    p = phyto.p
    initial_biomass = sum(phyto.data.Bm[phyto.data.ac])

    rich = deepcopy(initial)
    rich.tracers.NH4.data .= 3 * p.KsatNH4
    rich.tracers.NO3.data .= 0
    rich.tracers.PO4.data .= 3 * p.KsatPO4
    sim_run!(forced_example(rich, 200.0))
    reference = rich.individuals.phytos.sp1.data
    @test all(>(0), reference.PS[reference.ac])
    @test reference.BS[reference.ac] ≈ (3/4) .* reference.PS[reference.ac]
    @test sum(reference.Bm[reference.ac]) > initial_biomass

    @testset "$name" for (name, nh4, no3, po4, growth_fraction, nh4_share) in (
        ("NH4 only", p.KsatNH4, 0, 3 * p.KsatPO4, 1/2, 1.0),
        ("NO3 only", 0, p.KsatNO3, 3 * p.KsatPO4, 1/2, 0.0),
        ("Combined nitrogen", p.KsatNH4/2, p.KsatNO3/2, 3 * p.KsatPO4, 2/3, 1/2),
        ("Phosphate limitation", 3 * p.KsatNH4, 0, p.KsatPO4, 1/2, 1.0),
        ("Both at half saturation", p.KsatNH4, 0, p.KsatPO4, 1/2, 1.0),
        ("No inorganic nitrogen", 0, 0, 3 * p.KsatPO4, 0, 0.0),
        ("No phosphate", 3 * p.KsatNH4, 0, 0, 0, 1.0),
    )
        model = deepcopy(initial)
        model.tracers.NH4.data .= nh4
        model.tracers.NO3.data .= no3
        model.tracers.PO4.data .= po4
        inventory = carbon_mode_nutrient_inventory(model)
        carbon = example_inventory(model)
        sim = forced_example(model, 200.0)
        sim_run!(sim)
        data = model.individuals.phytos.sp1.data

        # Nutrients limit biosynthesis while the photosynthetic rate is unchanged.
        @test data.PS[data.ac] ≈ reference.PS[reference.ac]
        @test data.BS[data.ac] ≈ growth_fraction .* reference.PS[reference.ac]
        @test all(iszero, data.RP[data.ac])
        @test all(iszero, data.Bd[data.ac])
        @test data.Bm[data.ac] ≈ phyto.data.Bm[phyto.data.ac] .+ sim.ΔT .* data.BS[data.ac] rtol=1e-10

        # The tracer tendencies isolate particle exchange from background reactions.
        fixed = sim.ΔT * sum(data.BS[data.ac])
        flux = model.timestepper.plk
        @test tot_mass(flux.NH4.data, model.grid) ≈ -p.R_NC * fixed * nh4_share
        @test tot_mass(flux.NO3.data, model.grid) ≈ -p.R_NC * fixed * (1 - nh4_share)
        @test tot_mass(flux.PO4.data, model.grid) ≈ -p.R_PC * fixed
        @test carbon_mode_nutrient_inventory(model) ≈ inventory rtol=1e-10
        @test example_inventory(model) ≈ carbon rtol=1e-10
    end
end

@testset "CarbonMode repair with scarce nutrients" begin
    model = physiology_example(CarbonMode(); FT=Float64)
    data = model.individuals.phytos.sp1.data
    first, second = findall(data.ac)[1:2]
    data.Bd[first] = 0.2 * data.Bm[first]
    data.Bd[second] = 1e-6 * data.Bm[second]
    damage = copy(data.Bd)
    carbon = example_inventory(model)
    model.tracers.NH4.data .= 0
    model.tracers.NO3.data .= 0
    sim = forced_example(model, 200.0)
    sim_run!(sim)

    @test all(>(0), data.PS[data.ac])
    @test all(iszero, data.BS[data.ac])
    @test data.RP[first] ≈ data.PS[first]
    @test data.RP[second] * sim.ΔT ≈ damage[second]
    @test 0 <= data.Bd[first] < damage[first]
    @test data.Bd[second] ≈ 0 atol=eps(damage[second])
    @test all(iszero, data.RP[findall(data.ac)[3:end]])
    @test example_inventory(model) ≈ carbon rtol=1e-10
end

@testset "CarbonMode organic losses" begin
    @testset "$event" for (event, fractions) in ((:graz, (0.2, 0.4, 0.6)), (:mort, (0.3, 0.5, 0.7)))
        model = physiology_example(CarbonMode(); FT=Float64)
        phyto = model.individuals.phytos.sp1
        keys = Tuple(Symbol(event, :Frac, element) for element in (:C, :N, :P))
        phyto.p = merge(phyto.p, NamedTuple{keys}(fractions))
        data = phyto.data
        i = findfirst(data.ac)
        biomass = data.Bm[i]
        before = carbon_mode_nutrient_inventory(model)
        getproperty(data, event)[i] = 1
        sim_run!(forced_example(model, 0.0))

        @test !data.ac[i]
        flux = model.timestepper.plk
        for (dissolved, particulate, ratio, fraction) in
            zip((:DOC, :DON, :DOP), (:POC, :PON, :POP), (1, phyto.p.R_NC, phyto.p.R_PC), fractions)
            @test tot_mass(flux[dissolved].data, model.grid) ≈ biomass * ratio * fraction
            @test tot_mass(flux[particulate].data, model.grid) ≈ biomass * ratio * (1 - fraction)
        end
        @test carbon_mode_nutrient_inventory(model) ≈ before rtol=1e-10
    end
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
    don = model.tracers.DON
    set_bc!(don, CPU(); pos=:west, bc_c=1e-3, bc_v=1.0)
    set_bc!(don, CPU(); pos=:east, bc_c=fill(2e-4, 3, 4), bc_v=1.0)
    set_bc!(don, CPU(); pos=:south, bc_c=cat(fill(3e-4, 2, 4), fill(6e-4, 2, 4); dims=3), bc_v=1.0)
    set_bc!(don, CPU(); pos=:north, bc_c=cat(fill(1e-4, 2, 4), fill(2e-4, 2, 4); dims=3), bc_v=1.0)
    initial_DON = tot_mass(model.tracers.DON.data, model.grid)
    initial_C = example_inventory(model)
    sim_run!(forced_example(model, 0.0; iterations=2))
    # Signed inflow × face area × duration. Both time slices must be applied once.
    expected = (1e-3 - 2e-4) * (6 * 8) * 20 +
               ((3e-4 - 1e-4) + (6e-4 - 2e-4)) * (2 * 8) * 10
    @test tot_mass(model.tracers.DON.data, model.grid) - initial_DON ≈ expected rtol=1e-6
    @test example_inventory(model) ≈ initial_C rtol=5e-6
end
