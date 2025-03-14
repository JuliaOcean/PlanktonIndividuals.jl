function construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x   = zeros(FT, maxN), y   = zeros(FT, maxN), z   = zeros(FT, maxN),
                          xi  = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          Sz  = zeros(FT, maxN),
                          Bm  = zeros(FT, maxN), Bd  = zeros(FT, maxN), Chl = zeros(FT, maxN),
                          gen = zeros(FT, maxN), age = zeros(FT, maxN),
                          ac  = zeros(FT, maxN), idx = zeros(Int,maxN),
                          PS  = zeros(FT, maxN), BS  = zeros(FT, maxN), RS  = zeros(FT, maxN),
                          TD  = zeros(FT, maxN), RP  = zeros(FT, maxN),
                          graz= zeros(FT, maxN), mort= zeros(FT, maxN), dvid= zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :Cquota, :mean, :var, :α, :Φ, :Topt, :Tmax, :Ea, :PCmax, :Chl2C, :respir, :f_T2B,
                 :grz_P, :dvid_type, :dvid_P, :dvid_reg, :dvid_reg2, :mort_P, :mort_reg, :grazFracC, :mortFracC,
                 :thermal, :is_bact)

    pkeys = collect(keys(params))
    tmp = zeros(length(param_names))
    for i in 1:length(param_names)
        if length(findall(x->x==string(param_names[i]),pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(FT.(tmp))
    return plankton(data, p)
end

function generate_plankton!(plank, N::Int, g::AbstractGrid, arch::Architecture)
    mean = plank.p.mean
    var = plank.p.var
    Cquota = plank.p.Cquota
    Nsuper = plank.p.Nsuper
    Chl2C = plank.p.Chl2C

    plank.data.ac[1:N]  .= 1.0f0                                                                             # activity
    plank.data.gen[1:N] .= 1.0f0                                                                             # generation
    plank.data.age[1:N] .= 0.0f0                                                                             # age

    randn!(rng_type(arch), plank.data.Sz)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)

    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                         # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                         # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                         # z, unit: grid spacing, starting from 0
    plank.data.Sz  .= max.(1.0f0, plank.data.Sz .* var .+ mean) .* plank.data.ac                     # size
    plank.data.Bm  .= Cquota .* plank.data.Sz .* Nsuper                                              # Bm
    plank.data.Chl .= plank.data.Bm .* Chl2C                                                         # Chl

    mask_individuals!(plank.data, g, N, arch)
end
