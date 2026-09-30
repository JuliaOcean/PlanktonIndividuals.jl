using PlanktonKernels.Biogeochemistry: bgc_params_default, update_bgc_params
using PlanktonIndividuals.Individuals: phyt_params_default, colony_params_default, abiotic_params_default,
                                      update_phyt_params, update_colony_params, update_abiotic_params

function test_parameters()
    for f in (update_phyt_params, update_colony_params, update_abiotic_params)
        @test parentmodule(f) === PlanktonIndividuals.Individuals
    end
    tmp_bgc = Dict("κh" => 1f-5, "κv" => 1f-5)
    param_bgc = update_bgc_params(tmp_bgc, Float32)
    @test param_bgc["κh"] == 1f-5
    @test param_bgc["κv"] == 1f-5

    tmp_phyt_carbon = Dict("PCmax" => [1f-5])
    param_phyt_carbon = update_phyt_params(tmp_phyt_carbon, Float32; N = 1, mode = CarbonMode())
    @test param_phyt_carbon["PCmax"] == [1f-5]

    tmp_phyt_quota = Dict("PCmax" => [1f-5, 1f-5, 1f-5])
    param_phyt_quota = update_phyt_params(tmp_phyt_quota, Float32; N = 3, mode = QuotaMode())
    @test param_phyt_quota["PCmax"] == [1f-5, 1f-5, 1f-5]

    return nothing
end

function test_mode_parameter_dispatch()
    mode_modules = ((MacroMolecularMode(), PlanktonIndividuals.Individuals.MacroMolecular),
                    (IronEnergyMode(), PlanktonIndividuals.Individuals.IronEnergy),
                    (QuotaMode(), PlanktonIndividuals.Individuals.Quota),
                    (CarbonMode(), PlanktonIndividuals.Individuals.Carbon),
                    (ProteinMode(), PlanktonIndividuals.Individuals.Protein))
    for (mode, owner) in mode_modules
        @test which(phyt_params_default, (Int64, typeof(mode))).module === PlanktonIndividuals.Individuals
        one = phyt_params_default(1, mode)
        many = phyt_params_default(3, mode)
        @test keys(one) == keys(many)
        for name in keys(one)
            @test many[name] == fill(one[name][1], 3)
        end
        updated = update_phyt_params(Dict("Nsuper" => [2,3,4]), Float32; N=3, mode)
        @test updated["Nsuper"] == Float32[2,3,4]
        @test phyt_params_default(3,mode)["Nsuper"] == many["Nsuper"]
        @test_throws ArgumentError update_phyt_params(Dict("unknown"=>[1]),Float32;mode)
    end
    @test which(colony_params_default, (Int64,Vector{Int64},IronEnergyMode)).module === PlanktonIndividuals.Individuals
    colonies = update_colony_params([Dict("Nsuper"=>[2,3]),Dict("Nsuper"=>[4,5,6])],Float32;
                                    Ncl=2,Nsp=[2,3],mode=IronEnergyMode())
    @test colonies[1]["Nsuper"] == Float32[2,3]
    @test colonies[2]["Nsuper"] == Float32[4,5,6]
    @test which(abiotic_params_default, (Int64,)).module === PlanktonIndividuals.Individuals
    @test update_abiotic_params(Dict("Nsuper"=>[2,3]),Float32;N=2)["Nsuper"] == Float32[2,3]
    return nothing
end

@testset "Parameters" begin
    @testset "Parameter overrides" begin
        test_parameters()
    end
    @testset "Mode-dispatched defaults" begin
        test_mode_parameter_dispatch()
    end
end

function test_particle_diffusivity_parameters()
    inds = PlanktonIndividuals.Individuals
    for mode in (MacroMolecularMode(),IronEnergyMode(),QuotaMode(),CarbonMode(),ProteinMode())
        params = update_phyt_params(Dict("κhP"=>[0.0,0.02],"κvP"=>[0.0,0.03]),Float32;N=2,mode)
        @test phyt_params_default(2,mode)["κhP"] == [0.0,0.0]
        for sp in 1:2
            particle = inds.construct_plankton(CPU(),sp,params,8,Float32,mode)
            @test particle.p.κhP === Float32(params["κhP"][sp])
            @test particle.p.κvP === Float32(params["κvP"][sp])
        end
    end
    @test !haskey(bgc_params_default(),"κhP")
    @test !haskey(bgc_params_default(),"κvP")
    @test_throws ArgumentError update_bgc_params(Dict("κhP"=>0.1),Float32)
    params=update_abiotic_params(Dict("κhP"=>[0.02],"κvP"=>[0.03]),Float32)
    particle=inds.Abiotic.construct_abiotic_particle(CPU(),1,params,8,Float32)
    @test particle.p.κhP === 0.02f0
    @test particle.p.κvP === 0.03f0
    params=update_colony_params([Dict("κhP"=>[0.02],"κvP"=>[0.03])],Float32)
    particle=inds.ColonyIronEnergy.construct_plankton(CPU(),1,params[1],8,Float32)
    @test particle.p.κhP === 0.02f0
    @test particle.p.κvP === 0.03f0
    grid=RectilinearGrid(size=(4,4,4),x=(0,4),y=(0,4),z=(0,-4))
    model=PlanktonModel(CPU(),grid;mode=CarbonMode(),options=let opt = model_options(); opt.max_individuals = 32; opt end,
        phyto=phyto_setup(Dict("κhP"=>[0.0,0.01],"κvP"=>[0.0,0.01]),[8,8],2))
    positions=[copy(sp.data.x) for sp in model.individuals.phytos]
    sim=PlanktonSimulation(model;ΔT=1.0,iterations=1)
    update!(sim)
    @test model.individuals.phytos.sp1.data.x[1:8] == positions[1][1:8]
    @test model.individuals.phytos.sp2.data.x[1:8] != positions[2][1:8]
    return nothing
end

@testset "Particle diffusivities" begin
    test_particle_diffusivity_parameters()
end
