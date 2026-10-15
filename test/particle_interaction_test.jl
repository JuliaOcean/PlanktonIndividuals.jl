using PlanktonIndividuals, PlanktonKernels

@testset "Abiotic attachment across physiology modes" begin
    abiotic_module = PlanktonIndividuals.Individuals.Abiotic
    for mode in (CarbonMode(), QuotaMode(), MacroMolecularMode(), IronEnergyMode(), ProteinMode()),
        FT in (Float32, Float64)
        @testset "$(typeof(mode)), $FT" begin
            grid = RectilinearGrid(size=(1, 1, 1), x=(0, 1), y=(0, 1), z=(0, -1))
            model = PlanktonModel(CPU(), grid; mode, FT,
                                 options=model_options(max_individuals=64, max_candidates=2),
                                 phyto=phyto_setup(Dict("max_ptc"=>[3]), [2], 1),
                                 abiotic=abiotic_setup(Dict("Rd"=>[1.0]), [32], 1,
                                                       Palat([(:sp1, :sa1)], [])))
            phyto = model.individuals.phytos.sp1
            particles = model.individuals.abiotics.sa1
            @test phyto.p.max_ptc == 3
            @test PlanktonIndividuals.Individuals.phyt_params_default(1, mode)["max_ptc"] == [25]
            # Force a real encounter; sparse random examples rarely reach merge.
            for data in (phyto.data, particles.data)
                data.x[data.ac] .= 0.5
                data.y[data.ac] .= 0.5
                data.z[data.ac] .= 0.5
                data.xi .= grid.Hx + 1
                data.yi .= grid.Hy + 1
                data.zi .= grid.Hz + 1
            end
            intac = model.timestepper.intac
            abiotic_module.calc_interaction!(intac, phyto.data, particles.data,
                                             model.timestepper.rnd, model.grid,
                                             particles.p, 2, model.arch)
            @test all(>(0), intac[:, 1:32])
            @test all(iszero, intac[:, 33:end])
            before = count(particles.data.ac) + sum(phyto.data.ptc)
            abiotic_module.merge_particle!(intac, phyto.data, phyto.p, particles.data, 2, model.arch)
            @test phyto.data.ptc[1:2] == FT[3, 3]
            @test all(iszero, phyto.data.ptc[3:end])
            @test count(particles.data.ac) == 26
            @test count(particles.data.ac) + sum(phyto.data.ptc) == before
            # Full cells reject further attachment, and inactive particles cannot
            # be counted a second time even if their candidates are still present.
            saved_activity = copy(particles.data.ac)
            abiotic_module.merge_particle!(intac, phyto.data, phyto.p, particles.data, 2, model.arch)
            @test particles.data.ac == saved_activity
            @test phyto.data.ptc[1:2] == FT[3, 3]
            abiotic_module.reset_interaction!(intac, model.arch)
            @test all(iszero, intac)
            phyto.data.ptc .= 0
            abiotic_module.merge_particle!(intac, phyto.data, phyto.p, particles.data, 2, model.arch)
            @test all(iszero, phyto.data.ptc)
            @test particles.data.ac == saved_activity
            # Zero capacity is a supported way to disable attachment.
            intac[:, :] .= 1
            disabled = merge(phyto.p, (max_ptc=zero(FT),))
            abiotic_module.merge_particle!(intac, phyto.data, disabled, particles.data, 2, model.arch)
            @test all(iszero, phyto.data.ptc)
            @test particles.data.ac == saved_activity
        end
    end
end
