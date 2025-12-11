"""
    default_PARF(grid, ΔT, iterations)
Generate default hourly surface PAR.
"""
function default_PARF(grid, ΔT, iterations)
    PAR = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3871666666666666, 87.10258333333333, 475.78150000000016, 929.2737916666669,
           1232.3633333333337, 1638.918916666667, 1823.7921666666664, 1906.2769583333336, 1776.0280416666667,
           1678.5026249999999, 1410.216666666667, 815.4129583333336, 525.104, 135.993, 2.9493750000000003, 0.0, 0.0, 0.0,]
    total_days = Int((ΔT * iterations) ÷ 86400 + 1)
    PAR_day = zeros(grid.Nx, grid.Ny, 24)
    for i in 1:24
        PAR_day[:,:,i] .= PAR[i]
    end
    PAR_domain = repeat(PAR_day, outer = (1,1,total_days))
    return PAR_domain
end

"""
    default_temperature(grid, ΔT, iterations)
Generate default hourly temperature.
"""
function default_temperature(grid, ΔT, iterations)
    temp = [26.646446609406727, 26.56698729810778, 26.517037086855467, 26.5, 26.517037086855467, 26.56698729810778,
            26.646446609406727, 26.75, 26.87059047744874, 27.0, 27.12940952255126, 27.25, 27.353553390593273,
            27.43301270189222, 27.482962913144533, 27.5, 27.482962913144533, 27.43301270189222, 27.353553390593273,
            27.25, 27.12940952255126, 27.0, 26.87059047744874, 26.75]
    total_days = Int((ΔT * iterations) ÷ 86400 + 1)
    temp_day = zeros(grid.Nx, grid.Ny, grid.Nz, 24)
    for i in 1:24
        temp_day[:,:,end,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in grid.Nz-1:-1:1
        temp_day[:,:,j,:] .= temp_day[:,:,j+1,:] .+ 4.0e-4 * grid.zC[j]
    end
    temp_domain = repeat(temp_day, outer = (1,1,1,total_days))
    return temp_domain
end

"""
    bgc_params_default(FT::DataType)
Generate default biogeochemical parameter values 
"""
function bgc_params_default(FT)
    params = Dict{String, FT}(
        "kw"           => 0.046,              # PAR attenuation (/m)
        "kc"           => 0.04,               # PAR attenuation (m²/mgChl)
        "kDOC"         => 1/30/86400,         # Remineralization rate for DOC, turn over time: a month (per second)
        "Nit"          => 1/30/86400,         # Nitrification rate for NH4
        "kDON"         => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kDOP"         => 1/30/86400,         # Remineralization rate for DOP, turn over time: a month (per second)
        "kPOC"         => 1/30/86400,         # Remineralization rate for POC, turn over time: a month (per second)
        "kPON"         => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "kPOP"         => 1/30/86400,         # Remineralization rate for POP, turn over time: a month (per second)
        "κh"           => 0.0e-6,             # Horizontal diffusion
        "κv"           => 0.0e-6,             # Vertical diffusion
        "κhP"          => 0.0e-6,             # Horizontal diffusion for individuals
        "κvP"          => 0.0e-6,             # Vertical diffusion for individuals
        "shared_graz"  => 1.0,                # Using shared grazing or not - default is yes - 1.0
        "w_sink_org"   => 2.31e-8,            # Settling velocity of organic matter (m/s)
        "w_sink_inorg" => 2.31e-5,            # Settling velocity of inorganic matter (m/s)
        "lambda_POC"   => 5.8e-8,             # Scavenging rate of iron by POC (m³/mmol/second)
        "lambda_min"   => 3.5e-10,            # Minimum scavenging rate of iron (per second)
        "kdiss"        => 4.4e-7 ,           # Dissolution rate of inorganic particulate iron  (per second)
        "lambda_dust"  => 1.7e-3,             # Scavenging rate of iron by dust (m³/kg/second)
        "lambda_Fe"    => 1.7,                # Coagulation rate of dissolved iron (m³/mmol/second)
        "ligand"       => 6.0e-4,             # Ligand concentration (mmol/m³/second)
        "DFeFrac"      => 0.01,               # Free DFe Fraction
        )
    return params
end

"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
#=
CH includes cabohydrate and lipids
┌─────────┬────────────────┬────────────────────┬───────────────────────────┬──────────────────────────┐
│         │ Micromonas sp. │ Ostreococcus tauri │ Thalassiosira weissflogii │ Thalassiosira pseudonana │
├─────────┼────────────────┼────────────────────┼───────────────────────────┼──────────────────────────┤
│ PRO:DNA │          34.83 │              22.75 │                     67.99 │                   135.36 │
│ RNA:DNA │           1.73 │               1.71 │                      3.70 │                     6.49 │
│ CHL:DNA │           3.45 │               2.16 │                      9.83 │                    18.05 │
│  CH:DNA │          23.32 │              19.01 │                     79.26 │                   121.48 │
│   k_pro │       7.06e-05 │           1.50e-04 │                  7.64e-05 │                 1.17e-04 │
│   k_dna │       1.62e-07 │           5.44e-07 │                  1.16e-07 │                 9.61e-08 │
│   k_rna │       3.24e-07 │           6.71e-07 │                  6.60e-07 │                 6.60e-07 │
│k_pro_sat│       4.50e-13 │           4.82e-14 │                  3.21e-11 │                 2.73e-12 │
│k_dna_sat│       8.29e-16 │           2.07e-16 │                  1.45e-12 │                 6.01e-14 │
│k_rna_sat│       1.00e-12 │           2.40e-14 │                  1.76e-10 │                 1.74e-11 │
└─────────┴────────────────┴────────────────────┴───────────────────────────┴──────────────────────────┘
=#
function phyt_params_default(N::Int64, mode::MacroMolecularMode)
    params=Dict(
        "Nsuper"   => [1],       # Number of phyto cells each super individual represents
        "C_DNA"    => [1.8e-13], # DNA C quota of phyto cells (mmolC/cell)
        "var"      => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "RNA2DNA"  => [1.73],    # Initial RNA:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "PRO2DNA"  => [34.83],   # Initial protein:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "CH2DNA"   => [23.3],    # Initial carbohydrate+lipid:DNA ratio in phytoplankton (mmol C/mmolC) from Micromonas sp.
        "Chl2DNA"  => [3.45],    # Initial Chla:DNA ratio in phytoplankton (mmolC/mmolC) from Micromonas sp.
        "α"        => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"     => [27.0],    # Optimal temperature for growth (C)
        "Tmax"     => [30.0],    # Maximal temperature for growth (C)
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [6.2e-5],  # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6],  # Maximum P uptake rate (mmol P/mmol C/second)
        "KsatNH4"  => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0],     # Half-saturation coeff (mmol C/m³)
        "NSTmax"   => [0.12],    # Maximum N reserve in total N (mmol N/mmol N)
        "PSTmax"   => [0.70],    # Maximum P reserve in total P (mmol P/mmol P)
        "CHmax"    => [0.4],     # Maximum Carbohydrate in cell (mmol C/mmol C)
        "respir"   => [1.2e-6],  # Respiration rate(per second)
        "k_pro"    => [6.0e-5],  # Protein synthesis rate (mmol C/mmol C/second)
        "k_sat_pro"=> [4.5e-13], # Hafl saturation constent for protein synthesis (mmol C/cell)
        "k_dna"    => [1.6e-7],  # DNA synthesis rate (mmol C/mmol C/second)
        "k_sat_dna"=> [1.0e-15], # Hafl saturation constent for DNA synthesis (mmol C/cell)
        "k_rna"    => [3.0e-7],  # RNA synthesis rate (mmol C/mmol C/second)
        "k_sat_rna"=> [1.0e-12], # Hafl saturation constent for RNA synthesis (mmol C/cell)
        "Chl2N"    => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC_PRO" => [1/4.5],   # N:C ratio in protein (from Inomura et al 2020.)
        "R_NC_DNA" => [1/2.9],   # N:C ratio in DNA (from Inomura et al 2020.)
        "R_PC_DNA" => [1/11.1],  # P:C ratio in DNA
        "R_NC_RNA" => [1/2.8],   # N:C ratio in RNA (from Inomura et al 2020.)
        "R_PC_RNA" => [1/10.7],  # P:C ratio in RNA
        "dvid_P"   => [1.0e-5],  # Division probability per second
        "grz_P"    => [0.0],     # Grazing probability per second
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::IronEnergyMode)
    params=Dict(
        "Nsuper"    => [1],       # Number of phyto cells each super individual represents (cells)
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0 (mmolC/cell)
        "Rad"       => [0.12],    # Radius (μm) of Prochlorococcus
        "mean"      => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"       => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint"  => [0.10],    # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [4.5e-2],  # Irradiance absorption coeff (mmolC m² second/mgChl /μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea"        => [5.3e4],   # Free energy
        "is_nr"     => [1.0],     # 1 for non-diazotroph, 0 for diazotroph
        "is_croc"   => [0.0],     # 1 for Crocosphaera-like N fixation pattern
        "is_tric"   => [0.0],     # 1 for Trichodesmium-like N fixation pattern
        "PCmax"     => [8.0e-8],  # Maximum light harvesting rate (mmolATP/mmolC/second) 
        "VNH4max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VNO3max"   => [6.9e-6],  # Maximum N uptake rate (mmolN/mmolC/second)
        "VPO4max"   => [1.2e-6],  # Maximum P uptake rate (mmolP/mmolC/second)
        "k_O2"      => [3.0e-3],  # O₂ permeability coefficient (m/s)
        "k_cf"      => [1.0e-5],  # Carbon fixation rate (per second)
        "k_rs"      => [1.5e-6],  # Maximum respiration rate (per second)
        "k_nr"      => [2.8e-6],  # Nitrate reduction rate (per second)
        "k_nf"      => [2.8e-6],  # N fixation rate (mmolN/mmolC/second)
        "k_mtb"     => [3.5e-5],  # Metabolic rate (per second)
        "e_cf"      => [2.6],     # Energy consumption rate of carbon fixation (mmolATP/mmolC)
        "e_rs"      => [5.0],     # Energy production rate of respiration (mmolATP/mmolC)
        "e_nf"      => [8.0],     # Energy consumption rate of N fixation (mmolATP/mmolN)
        "e_nr"      => [0.0],     # Energy consumption rate of NO3 reduction (mmolATP/mmolN)
        "k_Fe_ST2PS"=> [2.4e-5],  # Allocation rate of Fe from storage to PS (per second)
        "k_Fe_PS2ST"=> [1.2e-6],  # Allocation rate of Fe from PS to storage (per second)
        "k_Fe_ST2NR"=> [1.2e-5],  # Allocation rate of Fe from storage to NR (per second)
        "k_Fe_NR2ST"=> [1.2e-5],  # Allocation rate of Fe from NR to storage (per second)
        "k_Fe_ST2NF"=> [1.2e-5],  # Allocation rate of Fe from storage to NF (per second)
        "k_Fe_NF2ST"=> [1.2e-5],  # Allocation rate of Fe from NF to storage (per second)
        "KfePS"     => [3.0e-6],  # Haff-saturation coeff of iron quota for photosynthesis (mmolFe/mmolC)
        "KfeNR"     => [2.0e-6],  # Haff-saturation coeff of iron quota for NO3 reduction (mmolFe/mmolC)
        "KfeNF"     => [5.0e-6],  # Haff-saturation coeff of iron quota for N fixation (mmolFe/mmolC)
        "KsatNH4"   => [0.005],   # Half-saturation coeff (mmolN/m³)
        "KsatNO3"   => [0.010],   # Half-saturation coeff (mmolN/m³)
        "KsatPO4"   => [0.003],   # Half-saturation coeff (mmolP/m³)
        "KSAFe"     => [2.77e-19],# Surface-area specific iron uptake rate (m³/μm²/cell/second)
        "qNO3max"   => [0.25],    # Maximum NO3 quota in cell (mmolN/mmolC)
        "qNH4max"   => [0.25],    # Maximum NH4 quota in cell (mmolN/mmolC)
        "qPmax"     => [0.02],    # Maximum P quota in cell (mmolP/mmolC)
        "qFemax"    => [2.0e-5],  # Maximum Fe quota in cell (mmolFe/mmolC)
        "CHmax"     => [0.4],     # Maximum C quota in cell (mmolC/mmolC)
        "qOmax"     => [0.5],     # Maximum O₂ quota in cell (mmolO₂/mmolC)
        "Chl2N"     => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"      => [16/106],  # N:C ratio in cell biomass
        "R_PC"      => [1/106],   # N:C ratio in cell biomass
        "NF_clock"  => [21600.0], # the circadian clock for N fixation
        "grz_P"     => [0.0],     # Grazing probability per second
        "dvid_P"    => [1e-4],    # Probability of cell division per second.
        "dvid_type" => [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg"  => [2.5],     # Regulations of cell division (cell size)
        "dvid_reg2" => [12.0],    # Regulations of cell division (clock time)
        "mort_P"    => [5e-5],    # Probability of cell natural death per second
        "mort_reg"  => [0.5],     # Regulation of cell natural death
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP" => [0.7],     # Fraction goes into dissolved organic pool
        "grazFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP" => [0.5],     # Fraction goes into dissolved organic pool
        "mortFracFe"=> [0.1],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::QuotaMode)
    params=Dict(
        "Nsuper"   => [1],       # Number of phyto cells each super individual represents
        "Cquota"   => [1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"     => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"      => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2Cint" => [0.10],    # Initial Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"        => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"     => [27.0],    # Optimal temperature for growth (C)
        "Tmax"     => [30.0],    # Maximal temperature for growth (C)
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [4.2e-5],  # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6],  # Maximum P uptake rate (mmol P/mmol C/second)
        "KsatNH4"  => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0],     # Half-saturation coeff (mmol C/m³)
        "Nqmax"    => [0.25],    # Maximum N quota in cell (mmol N/mmol C)
        "Pqmax"    => [0.02],    # Maximum P quota in cell (mmol P/mmol C)
        "Cqmax"    => [0.4],     # Maximum C quota in cell (mmol C/mmol C)
        "k_mtb"    => [3.5e-5],  # Metabolic rate (per second)
        "respir"   => [1.2e-6],  # Respiration rate(per second)
        "Chl2N"    => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"     => [16/106],  # N:C ratio in cell biomass
        "R_PC"     => [1/106],   # N:C ratio in cell biomass
        "grz_P"    => [0.0],     # Grazing probability per second
        "dvid_P"   => [1e-4],    # Probability of cell division per second.
        "dvid_type"=> [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg" => [2.5],     # Regulations of cell division (cell size)
        "dvid_reg2"=> [12.0],    # Regulations of cell division (clock time)
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    phyt_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyt_params_default(N::Int64, mode::CarbonMode)
    params=Dict(
        "Nsuper"    => [1],       # Number of phyto cells each super individual represents
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"      => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"       => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2C"     => [0.10],    # Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"         => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea"        => [5.3e4],   # Free energy
        "PCmax"     => [4.2e-5],  # Maximum primary production rate (per second)
        "respir"    => [1.2e-6],  # Respiration rate(per second)
        "f_T2B"     => [2.7e-6],  # Thermal damage rate (1.0/K/s)
        "grz_P"     => [0.0],     # Grazing probability per second
        "dvid_P"    => [5e-5],    # Probability of cell division per second.
        "dvid_type" => [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_reg"  => [1.9],     # Regulations of cell division (sizer)
        "dvid_reg2" => [12.0],    # Regulations of cell division (sizer)
        "mort_P"    => [5e-5],    # Probability of cell natural death per second
        "mort_reg"  => [0.5],     # Regulation of cell natural death
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
        "thermal"   => [1.0],     # thermal damage, 1 for on, 0 for off
        "is_bact"   => [0.0],     # is this specis bacteria or not, 1 for bacteria, 0 for phytoplankton
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    abiotic_params_default(N::Int64)
Generate default abiotic particle parameter values based on species number `N`.
"""
function abiotic_params_default(N::Int64)
    params=Dict(
        "Nsuper"    => [1],       # Number of abiotic particles each super individual represents
        "Rd"        => [1.0e-4],  # Distance between abiotic particle and phytoplankton cell (m)
        "release_P" => [1.0e-6],  # Probability of particle release per second
        "sz_min"    => [1.0e-3],  # Minimal size of a abiotic particle (mmolFe/particle)
        "Ktr"       => [0.0e-5],  # Rate of particle forming from tracer
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

function generate_n_species_params(N, params)
    p = []
    for key in keys(params)
        pa = (key, fill(params[key][1],N))
        push!(p, pa)
    end
    return Dict(p)
end