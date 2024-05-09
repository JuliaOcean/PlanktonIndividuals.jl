function construct_plankton(arch::Architecture, sp::Int, params::Dict, maxN::Int, FT::DataType)
    rawdata = StructArray(x    = zeros(FT, maxN), y    = zeros(FT, maxN), z    = zeros(FT, maxN),
                          xi   = zeros(Int,maxN), yi   = zeros(Int,maxN), zi   = zeros(Int,maxN),
                          CH   = zeros(FT, maxN), NST  = zeros(FT, maxN), PST  = zeros(FT, maxN),
                          PRO  = zeros(FT, maxN), DNA  = zeros(FT, maxN), RNA  = zeros(FT, maxN), 
                          Chl  = zeros(FT, maxN), gen  = zeros(FT, maxN), age  = zeros(FT, maxN), 
                          ac   = zeros(FT, maxN), idx  = zeros(FT, maxN),
                          PS   = zeros(FT, maxN), VDOC = zeros(FT, maxN), VNH4 = zeros(FT, maxN),
                          VNO3 = zeros(FT, maxN), VPO4 = zeros(FT, maxN), ρChl = zeros(FT, maxN),
                          resp = zeros(FT, maxN), S_PRO= zeros(FT, maxN), S_DNA= zeros(FT, maxN),
                          S_RNA= zeros(FT, maxN), exu  = zeros(FT, maxN),
                          graz = zeros(FT, maxN), mort = zeros(FT, maxN), dvid = zeros(FT, maxN)
                          ) 
    data = replace_storage(array_type(arch), rawdata)

    param_names=(:Nsuper, :C_DNA, :var, :PRO2DNA, :RNA2DNA, :CH2DNA, :Chl2DNA,
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
    p = NamedTuple{param_names}(FT.(tmp))
    return plankton(data, p)
end

function generate_plankton!(plank, N::Int, g::AbstractGrid, arch::Architecture)
    var = plank.p.var
    C_DNA = plank.p.C_DNA
    Nsuper = plank.p.Nsuper
    RNA2DNA = plank.p.RNA2DNA
    PRO2DNA = plank.p.PRO2DNA
    CH2DNA  = plank.p.CH2DNA
    Chl2DNA = plank.p.Chl2DNA

    plank.data.ac[1:N]  .= 1.0f0                    # activity
    plank.data.gen[1:N] .= 1.0f0                    # generation
    plank.data.age[1:N] .= 0.0f0                    # age

    rand!(rng_type(arch), plank.data.DNA)
    rand!(rng_type(arch), plank.data.RNA)
    rand!(rng_type(arch), plank.data.PRO)
    rand!(rng_type(arch), plank.data.x)
    rand!(rng_type(arch), plank.data.y)
    rand!(rng_type(arch), plank.data.z)
    rand!(rng_type(arch), plank.data.CH)
    rand!(rng_type(arch), plank.data.NST)
    rand!(rng_type(arch), plank.data.PST)

    plank.data.DNA .= plank.data.DNA .* var .+ 1.0f0                  # range: (1.0,1.0+var)
    plank.data.RNA .= plank.data.RNA .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PRO .= plank.data.PRO .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.CH  .= plank.data.CH  .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.NST .= plank.data.NST .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)
    plank.data.PST .= plank.data.PST .* var .* 2.0f0 .+ 1.0f0 .- var  # range: (1.0-var,1.0+var)


    plank.data.x   .=(plank.data.x .* g.Nx) .* plank.data.ac                                             # x, unit: grid spacing, starting from 0
    plank.data.y   .=(plank.data.y .* g.Ny) .* plank.data.ac                                             # y, unit: grid spacing, starting from 0
    plank.data.z   .=(plank.data.z .* g.Nz) .* plank.data.ac                                             # z, unit: grid spacing, starting from 0
    plank.data.DNA .= plank.data.DNA .* C_DNA .* Nsuper .* plank.data.ac                                 # DNA mmolC/individual
    plank.data.RNA .= plank.data.RNA .* C_DNA .* Nsuper .* plank.data.ac .* RNA2DNA                      # RNA mmolC/individual
    plank.data.PRO .= plank.data.PRO .* C_DNA .* Nsuper .* plank.data.ac .* PRO2DNA                      # PRO mmolC/individual
    plank.data.CH  .= plank.data.CH  .* plank.data.DNA .* CH2DNA                                         # CH  mmolC/individual
    plank.data.NST .= plank.data.NST .* plank.data.CH .* 16.0f0 ./ 106.0f0                               # NST mmolN/individual
    plank.data.PST .= plank.data.PST .* plank.data.CH ./ 106.0f0                                         # PST mmolP/individual
    plank.data.Chl .= plank.data.DNA .* Chl2DNA * 893.49f0 / 55.0f0                                      # Chl mgChl/individual

    mask_individuals!(plank.data, g, N, arch)
end

@inline function total_C_biomass(PRO, DNA, RNA, CH, Chl)
    C_tot = PRO + DNA + RNA + CH + Chl / 893.49f0 * 55.0f0 
    return C_tot
end
@inline function total_N_biomass(PRO, DNA, RNA, NST, Chl, p)
    N_tot = PRO * p.R_NC_PRO + DNA * p.R_NC_DNA + RNA * p.R_NC_RNA + NST + Chl / 893.49f0 * 4.0f0
    return N_tot
end
@inline function total_P_biomass(DNA, RNA, PST, p)
    P_tot = DNA * p.R_PC_DNA + RNA * p.R_PC_RNA + PST
    return P_tot
end
