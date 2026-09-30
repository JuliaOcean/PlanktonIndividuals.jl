@testset "Model options" begin
    opt=model_options()
    @test (opt.max_individuals,opt.max_candidates,opt.kc,opt.kw,opt.shared_graz) == (8192,25,0.04,0.046,1.0)
    opt.max_individuals=32
    opt.max_candidates=7
    opt.kc=0.05
    opt.kw=0.06
    opt.shared_graz=0.0
    grid=RectilinearGrid(size=(2,2,2),x=(0,2),y=(0,2),z=(0,-2))
    model=PlanktonModel(CPU(),grid;options=opt,mode=CarbonMode(),phyto=phyto_setup(nothing,[8],1))
    @test model.options === opt
    @test model.options.kc == 0.05
    @test model.options.kw == 0.06
    @test length(model.individuals.phytos.sp1.data.x) == 32
    @test model.options.max_candidates == 7
    @test !haskey(model.bgc_params,"shared_graz")
    update!(PlanktonSimulation(model;ΔT=1.0,iterations=1))
    @test all(isfinite,model.tracers.DIC.data)
    other_opt=deepcopy(opt)
    other_opt.max_individuals=64
    other_opt.max_candidates=9
    overridden=PlanktonModel(CPU(),grid;options=other_opt,
                             mode=CarbonMode(),phyto=phyto_setup(nothing,[8],1))
    @test overridden.options.max_individuals == 64
    @test overridden.options.max_candidates == 9
    @test (opt.max_individuals,opt.max_candidates) == (32,7)
end
