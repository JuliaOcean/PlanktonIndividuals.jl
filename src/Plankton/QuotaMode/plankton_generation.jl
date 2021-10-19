function construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)
    rawdata = StructArray(x    = zeros(maxN), y    = zeros(maxN), z    = zeros(maxN),
                          xi   = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          iS   = zeros(maxN), Sz   = zeros(maxN), Bm   = zeros(maxN), 
                          Cq   = zeros(maxN), Nq   = zeros(maxN), Pq   = zeros(maxN), 
                          Chl  = zeros(maxN), gen  = zeros(maxN), age  = zeros(maxN), 
                          ac   = zeros(maxN), idx  = zeros(maxN),
                          PS   = zeros(maxN), VDOC = zeros(maxN), VNH4 = zeros(maxN),
                          VNO3 = zeros(maxN), VPO4 = zeros(maxN), ρChl = zeros(maxN),
                          resp = zeros(maxN), BS   = zeros(maxN), exu  = zeros(maxN),
                          graz = zeros(maxN), mort = zeros(maxN), dvid = zeros(maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :Cquota, :mean, :var, :Chl2Cint, :α, :Φ, :T⁺, :Ea,
                 :PCmax, :PC_b, :VDOCmax, :VDOC_b, :VNO3max, :VNH4max, :VN_b, :VPO4max, :VP_b,
                 :KsatDOC, :KsatNH4, :KsatNO3, :KsatPO4, :Cqmax, :Cqmin, :Nqmax, :Nqmin, :Pqmax, :Pqmin,
                 :Chl2N, :R_NC, :R_PC, :k_mtb, :k_mtb_b, :respir_a, :respir_b,
                 :grz_P, :dvid_type, :dvid_P, :dvid_stp, :dvid_reg, :dvid_stp2, :dvid_reg2,
                 :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP, :mortFracC, :mortFracN, :mortFracP)

    pkeys = collect(keys(params))
    tmp = zeros(length(param_names))
    for i in 1:length(param_names)
        if length(findall(x->x==string(param_names[i]),pkeys))==0
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(tmp)
    return plankton(data, p)
end

function generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture; mask = nothing)
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

    plank.data.ac[1:N]  .= 1.0                                                                             # activity
    plank.data.gen[1:N] .= 1.0                                                                             # generation
    plank.data.age[1:N] .= 0.0                                                                             # age

    randn!(rng_type(arch), plank.data.iS)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.Cq)
    rand!(rng_type(arch), plank.data.Nq)
    rand!(rng_type(arch), plank.data.Pq)

    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                         # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                         # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                         # z, unit: grid spacing, starting from 0
    plank.data.iS  .= max.(1.0, plank.data.iS .* var .+ mean) .* plank.data.ac                       # init_size
    plank.data.Sz  .= copy(plank.data.iS)                                                            # size
    plank.data.Bm  .= Cquota .* plank.data.Sz .* Nsuper                                              # Bm
    plank.data.Cq  .=(plank.data.Cq .* (cqmax - cqmin)  .+ cqmin) .* plank.data.Bm                   # Cq
    plank.data.Nq  .=(plank.data.Nq .* (nqmax - nqmin)  .+ nqmin) .* plank.data.Bm                   # Nq
    plank.data.Pq  .=(plank.data.Pq .* (pqmax - pqmin)  .+ pqmin) .* plank.data.Bm                   # Pq
    plank.data.Chl .= plank.data.Bm .* Chl2Cint                                                      # Chl

    if mask ≠ nothing
        if size(mask) == (g.Nx, g.Ny, g.Nz)
            mask_individuals!(plank.data, mask, N, arch)
        else
            throw(ArgumentError("nut_mask: grid mismatch, size(mask) must equal to (grid.Nx, grid.Ny, grid.Nz)."))
        end
    end
end