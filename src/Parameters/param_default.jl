"""
    default_PARF(grid, ΔT, iterations)
Generate default hourly surface PAR.
"""
function default_PARF(grid, ΔT, iterations)
    PAR = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3871666666666666, 87.10258333333333, 475.78150000000016, 929.2737916666669,
           1232.3633333333337, 1638.918916666667, 1823.7921666666664, 1906.2769583333336, 1776.0280416666667,
           1678.5026249999999, 1410.216666666667, 815.4129583333336, 525.104, 135.993, 2.9493750000000003, 0.0, 0.0, 0.0,]
    total_days = (ΔT * iterations) ÷ 86400 + 1
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
    total_days = (ΔT * iterations) ÷ 86400 + 1
    temp_day = zeros(grid.Nx, grid.Ny, grid.Nz, 24)
    for i in 1:24
        temp_day[:,:,end,i] .= temp[i]
    end
    # vertical temperature gradient
    for j in grid.Nz-1:-1:1
        temp_day[:,:,j,:] .= temp_day[:,:,j+1,:] .- (0.04*(grid.zC[j+1]-grid.zC[j]))
    end
    temp_domain = repeat(temp_day, outer = (1,1,1,total_days))
    return temp_domain
end

"""
    bgc_params_default()
Generate default biogeochemical parameter values 
"""
function bgc_params_default()
    params = Dict(
        "kw"       => 0.046,              # PAR attenuation (/m)
        "kc"       => 0.04,               # PAR attenuation (m²/mgChl)
        "kDOC"     => 1/30/86400,         # Remineralization rate for DOC, turn over time: a month (per second)
        "Nit"      => 1/30/86400,         # Nitrification rate for NH4
        "kDON"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kDOP"     => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kPOC"     => 1/30/86400,         # Remineralization rate for POC, turn over time: a month (per second)
        "kPON"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "kPOP"     => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "κh"       => 0.0e-6,             # Horizontal diffusion
        "κv"       => 0.0e-6,             # Vertical diffusion
        "κhP"      => 0.0e-6,             # Horizontal diffusion for individuals
        "κvP"      => 0.0e-6,             # Vertical diffusion for individuals
    )
    return params
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
        "T⁺"       => [305.0],   # Maximal temperature for growth
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [4.2e-5],  # Maximum primary production rate (per second)
        "VDOCmax"  => [0.0],     # Maximum DOC uptake rate (mmol C/mmol C/second)
        "VNH4max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VNO3max"  => [6.9e-6],  # Maximum N uptake rate (mmol N/mmol C/second)
        "VPO4max"  => [1.2e-6],  # Maximum P uptake rate (mmol P/mmol C/second)
        "PC_b"     => [0.6],     # Shape parameter for size
        "VDOC_b"   => [0.6],     # Shape parameter for size
        "VN_b"     => [0.6],     # Shape parameter for size
        "VP_b"     => [0.6],     # Shape parameter for size
        "KsatNH4"  => [0.005],   # Half-saturation coeff (mmol N/m³)
        "KsatNO3"  => [0.010],   # Half-saturation coeff (mmol N/m³)
        "KsatPO4"  => [0.003],   # Half-saturation coeff (mmol P/m³)
        "KsatDOC"  => [0.0],     # Half-saturation coeff (mmol C/m³)
        "Nqmax"    => [0.12],    # Maximum N quota in cell (mmol N/mmol C)
        "Nqmin"    => [0.05],    # Minimum N quota in cell (mmol N/mmol C)
        "Pqmax"    => [0.01],    # Maximum P quota in cell (mmol P/mmol C)
        "Pqmin"    => [0.004],   # Minimum P quota in cell (mmol P/mmol C)
        "Cqmax"    => [0.4],     # Maximum C quota in cell (mmol C/mmol C)
        "Cqmin"    => [0.1],     # Minimum C quota in cell (mmol C/mmol C)
        "k_mtb"    => [3.5e-5],  # Metabolic rate (per second)
        "k_mtb_b"  => [0.25],    # Metabolic rate
        "respir_a" => [1.2e-6],  # Respiration rate(per second)
        "respir_b" => [0.6],     # Shape parameter for size
        "Chl2N"    => [3.0],     # Maximum Chla:N ratio in phytoplankton
        "R_NC"     => [16/106],  # N:C ratio in cell biomass
        "R_PC"     => [1/106],   # N:C ratio in cell biomass
        "grz_P"    => [0.0],     # Grazing probability per second
        "dvid_P"   => [5e-5],    # Probability of cell division per second.
        "dvid_type"=> [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_stp" => [6.0],     # Steepness of sigmoidal function
        "dvid_reg" => [1.9],     # Regulations of cell division (cell size)
        "dvid_stp2"=> [2.0],     # Steepness of sigmoidal function
        "dvid_reg2"=> [12.0],    # Regulations of cell division (clock time)
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracN"=> [0.7],     # Fraction goes into dissolved organic pool
        "grazFracP"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracN"=> [0.5],     # Fraction goes into dissolved organic pool
        "mortFracP"=> [0.5],     # Fraction goes into dissolved organic pool
        "ther_mort"=> [0],       # thermal mortality, 1 for on, 0 for off
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

function phyt_params_default(N::Int64, mode::CarbonMode)
    params=Dict(
        "Nsuper"   => [1],       # Number of phyto cells each super individual represents
        "Cquota"   => [1.8e-11], # C quota of phyto cells at size = 1.0
        "mean"     => [1.2],     # Mean of the normal distribution of initial phyto individuals
        "var"      => [0.3],     # Variance of the normal distribution of initial phyto individuals
        "Chl2C"    => [0.10],    # Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"        => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"        => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "T⁺"       => [298.0],   # Maximal temperature for growth
        "Ea"       => [5.3e4],   # Free energy
        "PCmax"    => [4.2e-5],  # Maximum primary production rate (per second)
        "PC_b"     => [0.6],     # Shape parameter for size
        "respir_a" => [1.2e-6],  # Respiration rate(per second)
        "respir_b" => [0.6],     # Shape parameter for size
        "grz_P"    => [0.0],     # Grazing probability per second
        "dvid_P"   => [5e-5],    # Probability of cell division per second.
        "dvid_type"=> [1],       # The type of cell division, 1:sizer, 2:adder.
        "dvid_stp" => [6.0],     # Steepness of sigmoidal function
        "dvid_reg" => [1.9],     # Regulations of cell division (sizer)
        "dvid_stp2"=> [2.0],     # Steepness of sigmoidal function
        "dvid_reg2"=> [12.0],    # Regulations of cell division (sizer)
        "mort_P"   => [5e-5],    # Probability of cell natural death per second
        "mort_reg" => [0.5],     # Regulation of cell natural death
        "grazFracC"=> [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC"=> [0.5],     # Fraction goes into dissolved organic pool
        "ther_mort"=> [0],       # thermal mortality, 1 for on, 0 for off
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