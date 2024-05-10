function construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x    = zeros(FT, maxN), y    = zeros(FT, maxN), z    = zeros(FT, maxN),
                          xi   = zeros(Int,maxN), yi   = zeros(Int,maxN), zi   = zeros(Int,maxN),
                          iS   = zeros(FT, maxN), Sz   = zeros(FT, maxN), Bm   = zeros(FT, maxN), 
                          Cq   = zeros(FT, maxN), Nq   = zeros(FT, maxN), Pq   = zeros(FT, maxN), 
                          Chl  = zeros(FT, maxN), gen  = zeros(FT, maxN), age  = zeros(FT, maxN), 
                          ac   = zeros(FT, maxN), idx  = zeros(Int, maxN),
                          PS   = zeros(FT, maxN), VDOC = zeros(FT, maxN), VNH4 = zeros(FT, maxN),
                          VNO3 = zeros(FT, maxN), VPO4 = zeros(FT, maxN), ρChl = zeros(FT, maxN),
                          resp = zeros(FT, maxN), BS   = zeros(FT, maxN), exu  = zeros(FT, maxN),
                          graz = zeros(FT, maxN), mort = zeros(FT, maxN), dvid = zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :Cquota, :mean, :var, :Chl2Cint, :α, :Φ, :Topt, :Tmax, :Ea,
                 :PCmax, :VDOCmax, :VNO3max, :VNH4max, :VPO4max,
                 :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :Cqmax, :Cqmin, :Nqmax, :Nqmin, :Pqmax, :Pqmin,
                 :Chl2N, :R_NC, :R_PC, :k_mtb, :respir, :grz_P, :dvid_type, :dvid_P, :dvid_reg, :dvid_reg2,
                 :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :mortFracC, :mortFracN, :mortFracP, :ther_mort)

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
    cqmax = plank.p.Cqmax
    cqmin = plank.p.Cqmin
    nqmax = plank.p.Nqmax
    nqmin = plank.p.Nqmin
    pqmax = plank.p.Pqmax
    pqmin = plank.p.Pqmin
    Chl2Cint = plank.p.Chl2Cint

    plank.data.ac[1:N]  .= 1.0f0                                                      # activity
    plank.data.gen[1:N] .= 1.0f0                                                      # generation
    plank.data.age[1:N] .= 0.0f0                                                      # age

    randn!(rng_type(arch), plank.data.iS)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.Cq)
    rand!(rng_type(arch), plank.data.Nq)
    rand!(rng_type(arch), plank.data.Pq)

    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                          # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                          # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                          # z, unit: grid spacing, starting from 0
    plank.data.iS  .= max.(1.0f0, plank.data.iS .* var .+ mean) .* plank.data.ac      # init_size
    plank.data.Sz  .= copy(plank.data.iS)                                             # size
    plank.data.Bm  .= Cquota .* plank.data.Sz .* Nsuper                               # Bm
    plank.data.Cq  .=(plank.data.Cq .* (cqmax - cqmin)  .+ cqmin) .* plank.data.Bm    # Cq
    plank.data.Nq  .=(plank.data.Nq .* (nqmax - nqmin)  .+ nqmin) .* plank.data.Bm    # Nq
    plank.data.Pq  .=(plank.data.Pq .* (pqmax - pqmin)  .+ pqmin) .* plank.data.Bm    # Pq
    plank.data.Chl .= plank.data.Bm .* Chl2Cint                                       # Chl

    mask_individuals!(plank.data, g, N, arch)
end
