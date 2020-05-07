# Parameters
param_default=Dict(
    "P_Nind"   => 1000,               # number of phyto individuals of each species
    "P_Nsp"    => 1,                  # number of phyto species
    "P_Nsuper" => 1e10,               # number of phyto cells each super individual represents
    "P_Cquota" => [1.8e-11, 1.8e-10], # C quota of phyto cells at size = 1.0
    "useDOC"   => [0,0],              # C quota of phyto cells at size = 1.0
    "isDiaz"   => [0,0],              # C quota of phyto cells at size = 1.0
    "P_mean"   => 1.5,                # mean of the normal distribution of initial phyto individuals
    "P_var"    => 0.3,                # variance of the normal distribution of initial phyto individuals
    "Z_Nind"   => 0,                  # number of zoo individuals of each species
    "Z_Nsp"    => 1,                  # number of zoo species
    "Z_Nsuper" => 1e0,                # number of zoo cells each super individual represents
    "Z_Cquota" => 0.0,                # C quota of zoo cells at size = 1.0
    "Z_mean"   => 0.0,                # mean of the normal distribution of initial zoo individuals
    "Z_var"    => 0.0,                # variance of the normal distribution of initial zoo individuals
    "diag_inds"=> zeros(Int,17),      # diagnostic indices, refer to diagnostics.jl
    "diag_freq"=> 3600,               # frequency of diagnostics (second)
    "PCmax"    => [1.8, 1.8],         # Maximum primary production rate (per day)
    "PC_b"     => [0.6, 0.6],         # Shape parameter for size
    "Chl2N"    => 3.0,                # Maximum Chla:N ratio in phytoplankton
    "Chl2Cint" => 0.16,               # Initial Chla:C ratio in phytoplankton (gChl/molC)
    "R_NC"     => 16/106,             # N:C ratio in cell biomass
    "R_PC"     => 1/106,              # N:C ratio in cell biomass
    "α"        => 8.0e-7,             # Irradiance absorption coeff (m^2/gChl)
    "katten_w" => 0.046,              # PAR attenuation (/m)
    "katten_c" => 0.04,               # PAR attenuation (/mgChl/m^3/m)
    "inhibcoef"=> [0.0, 0.0],         # light inhibition
    "Grz_P"    => 4000,               # phyt.size/Grz_P is the probability to be grazed
    "dvid_type"=> [1,1],              # the type of cell division, 1:sizer, 2:adder.
    "dvid_size"=> 1.9,                # relative cell size a cell can start divide
    "dvid_stp" => 6.0,                # steepness of sigmoidal function
    "dvid_add" => 1.5,                # added cell size for adder-like cell division
    "death_age"=> 168.0,              # average life time of a cell (7*24 hours)
    "a_death"  => 0.05,               # shape parameter for natural death
    "b_death"  => 0.5,                # shape parameter for natural death
    "Tempref"  => 293.15,             # reference temperature in K
    "TempAe"   => -4000.0,            # Arrenhius equation
    "TempCoeff"=> 0.8,                # Arrenhius equation
    "VNH4max"  => [0.6, 0.6],         # Maximum N uptake rate (mmol N/mmol C/day)
    "VNO3max"  => [0.6, 0.6],         # Maximum N uptake rate (mmol N/mmol C/day)
    "VPmax"    => [0.1, 0.1],         # Maximum P uptake rate (mmol P/mmol C/day)
    "VDOCmax"  => [0.0, 1.0],         # Maximum DOC uptake rate (mmol C/mmol C/day)
    "Nqmax"    => 0.17,               # Maximum N quota in cell (mmol N/mmol C)
    "Nqmin"    => 0.05,               # Minimum N quota in cell (mmol N/mmol C)
    "Pqmax"    => 0.02,               # Maximum P quota in cell (mmol P/mmol C)
    "Pqmin"    => 0.004,              # Minimum P quota in cell (mmol P/mmol C)
    "Cqmax"    => 0.4,                # Maximum C quota in cell (mmol C/mmol C)
    "Cqmin"    => 0.1,                # Minimum C quota in cell (mmol C/mmol C)
    "KsatNH4"  => [0.005, 0.005],     # Half-saturation coeff
    "KsatNO3"  => [0.010, 0.010],     # Half-saturation coeff
    "KsatP"    => [0.003, 0.003],     # Half-saturation coeff
    "KsatDOC"  => [0.0, 0.3],         # Half-saturation coeff
    "VN_b"     => [0.6, 0.6],         # Shape parameter for size
    "VP_b"     => [0.6, 0.6],         # Shape parameter for size
    "VDOC_b"   => [0.6, 0.6],         # Shape parameter for size
    "a_β"      => 3.1,                # scale parameter for metabolic partitioning of biosynthesis
    "b_β"      => -3.8,               # shape parameter for metabolic pratitioning of biosynthesis
    "k_mtb"    => 1.0/86400,          # metabolic rate (per second)
    "b_k_mtb"  => 0.5,                # metabolic rate (per second)
    "respir_ex"=> 0.0,                # Extra cost of C for biosynthesis
    "respir_b" => 0.13,               # Shape parameter for size
    "grazFracC"=> 0.7,                # Fraction goes into dissolved organic pool
    "grazFracN"=> 0.7,                # Fraction goes into dissolved organic pool
    "grazFracP"=> 0.7,                # Fraction goes into dissolved organic pool
    "mortFracC"=> 0.5,                # Fraction goes into dissolved organic pool
    "mortFracN"=> 0.5,                # Fraction goes into dissolved organic pool
    "mortFracP"=> 0.5,                # Fraction goes into dissolved organic pool
    "slpyFracC"=> 0.4,                # Fraction goes into zooplankton
    "slpyFracN"=> 0.4,                # Fraction goes into zooplankton
    "slpyFracP"=> 0.4,                # Fraction goes into zooplankton
    "k_sink"   => 0.01/86400,         # Sink rates for agents (m/s)
    "kDOC"     => 1/30/86400,         # remineralization rate for DOC, turn over time: a month (per second)
    "Nit"      => 1/10/86400,         # nitrification rate for NH4
    "kDON"     => 1/30/86400,         # remineralization rate for DON, turn over time: a month (per second)
    "kDOP"     => 1/30/86400,         # remineralization rate for DON, turn over time: a month (per second)
    "kPOC"     => 1/30/86400,         # remineralization rate for POC, turn over time: a month (per second)
    "kPON"     => 1/30/86400,         # remineralization rate for PON, turn over time: a month (per second)
    "kPOP"     => 1/30/86400,         # remineralization rate for PON, turn over time: a month (per second)
    "κh"       => 0,                  # horizontal diffusion
    "κv"       => 0,                  # vertical diffusion
    "κhP"      => 0,                  # horizontal diffusion for individuals
    "κvP"      => 0,                  # vertical diffusion for individuals
    "v_zoo"    => 1.0e-5              # velocity of zooplankton, m/s
)

# Options
#                    GridChoice, Gridoff, VelChoice, Veloff
RunOption=RunOptions(false,      Dict(),  false,     Dict())

#                  nTime, ΔT,   params,        Zoo
RunParam=RunParams(10,    600,  param_default, false)
