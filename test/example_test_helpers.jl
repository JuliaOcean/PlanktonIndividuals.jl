using Test, Random
using PlanktonIndividuals, PlanktonKernels

# Integrate concentrations over physical cells only, accumulating in Float64.
# Particle pools already contain amounts per super-individual (no volume factor).
function tot_mass(tracer, g)
    return sum(Float64(tracer[i+g.Hx, j+g.Hy, k+g.Hz]) *
               PlanktonKernels.Grids.volume(i+g.Hx, j+g.Hy, k+g.Hz, g)
               for i in 1:g.Nx, j in 1:g.Ny, k in 1:g.Nz)
end

function particle_phosphorus(phyto, mode)
    d, p = phyto.data, phyto.p
    if mode isa QuotaMode
        pool = d.Pq .+ d.Bm .* p.R_PC
    elseif mode isa IronEnergyMode
        pool = d.qP .+ d.Bm .* p.R_PC
    else # Molecular modes store P in nucleic acids and phosphorus reserves.
        pool = d.PST .+ d.DNA .* p.R_PC_DNA .+ d.RNA .* p.R_PC_RNA
    end
    return sum(Float64, pool[d.ac])
end

function example_inventory(model)
    g, tr = model.grid, model.tracers
    if model.mode isa CarbonMode
        return sum(tot_mass(tr[key].data, g) for key in (:DIC, :DOC, :POC)) +
               sum(sum(Float64, sp.data.Bm[sp.data.ac]) for sp in model.individuals.phytos)
    end
    return sum(tot_mass(tr[key].data, g) for key in (:PO4, :DOP, :POP)) +
           sum(particle_phosphorus(sp, model.mode) for sp in model.individuals.phytos)
end

function test_example_state(model; initial_counts=nothing)
    g = model.grid
    for (s, sp) in enumerate(model.individuals.phytos)
        d = sp.data
        @test length(d.ac) == model.options.max_individuals
        @test any(d.ac)
        if !isnothing(initial_counts)
            @test count(d.ac) == initial_counts[s]
            @test all(==(1), d.gen[d.ac])
            @test all(iszero, d.age[d.ac])
            @test all(iszero, d.Chl[.!d.ac])
            p = sp.p
            if model.mode isa CarbonMode || model.mode isa QuotaMode || model.mode isa IronEnergyMode
                @test d.Bm ≈ d.Sz .* (p.Cquota * p.Nsuper)
                ratio = model.mode isa CarbonMode ? p.Chl2C : p.Chl2Cint
                @test d.Chl ≈ d.Bm .* ratio
                if model.mode isa QuotaMode
                    @test all(d.Cq .<= p.Cqmax .* d.Bm)
                    @test all(d.Nq .+ p.R_NC .* d.Bm .<= p.Nqmax .* (d.Bm .+ d.Cq))
                    @test all(d.Pq .+ p.R_PC .* d.Bm .<= p.Pqmax .* (d.Bm .+ d.Cq))
                elseif model.mode isa IronEnergyMode
                    total_Fe = d.qFe .+ d.qFePS .+ d.qFeNR .+ d.qFeNF
                    @test d.qFePS ≈ 0.4 .* total_Fe
                    @test d.qFeNR ≈ (0.3 * p.is_nr) .* total_Fe
                    @test d.qFeNF ≈ (0.3 * (p.is_croc + p.is_tric)) .* total_Fe
                end
            else
                @test all(x -> p.C_DNA * p.Nsuper <= x <= p.C_DNA * p.Nsuper * (1 + p.var), d.DNA[d.ac])
                if model.mode isa MacroMolecularMode
                    @test d.Chl ≈ d.DNA .* (p.Chl2DNA * 893.49f0 / 55.0f0)
                else
                    @test all(>(0), d.PRO_PS[d.ac])
                    @test all(>(0), d.Chl_lag)
                end
            end
        end
        # Coordinates are in grid units; these checks catch failed boundary handling.
        for (coord, extent) in ((d.x, g.Nx), (d.y, g.Ny), (d.z, g.Nz))
            @test all(isfinite, coord)
            @test all(x -> 0 <= x <= extent, coord[d.ac])
        end
        # Check active physical pools; inactive slots are scratch space and may
        # contain undefined ratios. Their exclusion from coupling is checked by
        # the tracer and inventory assertions.
        for key in (:Bm, :Bd, :Cq, :Nq, :Pq, :CH, :NST, :PST, :DNA, :RNA, :PRO,
                    :PRO_RB, :PRO_MC, :PRO_MN, :PRO_TN, :PRO_TP, :PRO_TFe,
                    :PRO_RS, :PRO_PS, :PRO_OT, :qNH4, :qNO3, :qP, :qFe,
                    :qFePS, :qFeNR, :qFeNF, :Chl)
            hasproperty(d, key) || continue
            @test all(isfinite, getproperty(d, key)[d.ac])
            @test all(>=(0), getproperty(d, key)[d.ac])
        end
        @test all(>(0), d.Chl[d.ac])
    end
    for tr in model.tracers
        interior = @view tr.data[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz]
        @test all(isfinite, interior)
        @test all(>=(0), interior)
    end
end

function test_example_evolution!(sim)
    model = sim.model
    inventory = example_inventory(model)
    initial_DOC = copy(model.tracers.DOC.data)
    initial_positions = copy(model.individuals.phytos.sp1.data.z)
    t0, iter0 = model.t, model.iteration
    update!(sim)
    @test model.iteration == iter0 + sim.iterations
    @test model.t ≈ t0 + sim.ΔT * sim.iterations
    test_example_state(model)
    # Relative tolerance allows Float32 accumulation, but rejects percent-scale leaks.
    @test isapprox(example_inventory(model), inventory; rtol=5e-6, atol=0)
    @test model.tracers.DOC.data != initial_DOC
    if model.individuals.phytos.sp1.p.κvP > 0
        @test model.individuals.phytos.sp1.data.z != initial_positions
    end
end
