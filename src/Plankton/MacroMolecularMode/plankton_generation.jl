function construct_plankton(arch::Architecture, sp::Int64, params::Dict, maxN)
    rawdata = StructArray(x    = zeros(maxN), y    = zeros(maxN), z    = zeros(maxN),
                          xi   = zeros(Int,maxN), yi  = zeros(Int,maxN), zi  = zeros(Int,maxN),
                          CH   = zeros(maxN), NST  = zeros(maxN), PST  = zeros(maxN),
                          PRO  = zeros(maxN), DNA  = zeros(maxN), RNA  = zeros(maxN), 
                          Chl  = zeros(maxN), gen  = zeros(maxN), age  = zeros(maxN), 
                          ac   = zeros(maxN), idx  = zeros(maxN),
                          PS   = zeros(maxN), VDOC = zeros(maxN), VNH4 = zeros(maxN),
                          VNO3 = zeros(maxN), VPO4 = zeros(maxN), ρChl = zeros(maxN),
                          resp = zeros(maxN), S_PRO= zeros(maxN), S_DNA= zeros(maxN),
                          S_RNA= zeros(maxN), exu  = zeros(maxN),
                          graz = zeros(maxN), mort = zeros(maxN), dvid = zeros(maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :C_DNA, :mean, :var, :PRO2DNA, :RNA2DNA, :CH2DNA, :Chl2DNA,
                 :α, :Φ, :Topt, :Tmax, :Ea,
                 :PCmax, :VDOCmax, :VNO3max, :VNH4max, :VPO4max, :KsatDOC, :KsatNH4, :KsatNO3,
                 :KsatPO4, :CHmax, :CHmin, :NSTmax, :NSTmin, :PSTmax, :PSTmin,
                 :Chl2N, :R_NC_PRO, :R_NC_DNA, :R_NC_RNA, :R_PC_DNA, :R_PC_RNA, :respir,
                 :k_pro, :k_sat_pro, :k_rna, :k_sat_rna, :k_dna, :k_sat_dna, 
                 :dvid_P, :grz_P, :mort_P, :mort_reg, :grazFracC, :grazFracN, :grazFracP,
                 :mortFracC, :mortFracN, :mortFracP)

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

function generate_plankton!(plank, N::Int64, g::AbstractGrid, arch::Architecture)
    mean = plank.p.mean
    var = plank.p.var
    C_DNA = plank.p.C_DNA
    Nsuper = plank.p.Nsuper
    RNA2DNA = plank.p.RNA2DNA
    PRO2DNA = plank.p.PRO2DNA
    CH2DNA  = plank.p.CH2DNA
    Chl2DNA = plank.p.Chl2DNA

    plank.data.ac[1:N]  .= 1.0                                                                             # activity
    plank.data.gen[1:N] .= 1.0                                                                             # generation
    plank.data.age[1:N] .= 0.0                                                                             # age

    rand!(rng_type(arch), plank.data.DNA)
    rand!(rng_type(arch), plank.data.RNA)
    rand!(rng_type(arch), plank.data.PRO)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.CH)
    rand!(rng_type(arch), plank.data.NST)
    rand!(rng_type(arch), plank.data.PST)

    plank.data.DNA .= plank.data.DNA .* var .+ 1.0                # range: (1.0,1.0+var)
    plank.data.RNA .= plank.data.RNA .* var .* 2.0 .+ 1.0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO .= plank.data.PRO .* var .* 2.0 .+ 1.0 .- var  # range: (1.0-var,1.0+var)
    plank.data.CH  .= plank.data.CH  .* var .* 2.0 .+ 1.0 .- var  # range: (1.0-var,1.0+var)
    plank.data.NST .= plank.data.NST .* var .* 2.0 .+ 1.0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PST .= plank.data.PST .* var .* 2.0 .+ 1.0 .- var  # range: (1.0-var,1.0+var)


    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                         # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                         # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                         # z, unit: grid spacing, starting from 0
    plank.data.DNA .= plank.data.DNA .* C_DNA .* Nsuper .* plank.data.ac                             # DNA mmolC/individual
    plank.data.RNA .= plank.data.RNA .* C_DNA .* Nsuper .* plank.data.ac .* RNA2DNA                  # RNA mmolC/individual
    plank.data.PRO .= plank.data.PRO .* C_DNA .* Nsuper .* plank.data.ac .* PRO2DNA                  # PRO mmolC/individual
    plank.data.CH  .= plank.data.CH  .* plank.data.DNA .* CH2DNA                                     # CH  mmolC/individual
    plank.data.NST .= plank.data.NST .* plank.data.CH .* 16 ./ 106                                   # NST mmolN/individual
    plank.data.PST .= plank.data.PST .* plank.data.CH ./ 106                                         # PST mmolP/individual
    plank.data.Chl .= plank.data.DNA .* Chl2DNA * 893.49 / 55.0                                      # Chl mgChl/individual

    mask_individuals!(plank.data, g, N, arch)
end

@inline function total_C_biomass(PRO, DNA, RNA, CH, Chl)
    C_tot = PRO + DNA + RNA + CH + Chl / 893.49 * 55.0 
    return C_tot
end
@inline function total_N_biomass(PRO, DNA, RNA, NST, Chl, p)
    N_tot = PRO * p.R_NC_PRO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA + NST + Chl / 893.49 * 4.0
    return N_tot
end
@inline function total_P_biomass(DNA, RNA, PST, p)
    P_tot = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA + PST
    return P_tot
end
