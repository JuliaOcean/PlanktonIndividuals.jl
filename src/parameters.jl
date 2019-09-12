####################################### Parameters ########################################
PCmax = [1.8, 1.8] ./86400   # Maximum primary production rate (per second)
PC_b = [0.6, 0.6]            # Shape parameter for size
Chl2N = 3.0                  # Maximum Chla:N ratio in phytoplankton
R_NC  = 16/120               # N:C ratio in cell biomass, should be lower than 16:106
α = 8.0e-7                   # Irradiance absorption coeff (m^2/gChl)
katten_w = 0.046             # PAR attenuation (/m)
katten_c = 0.04              # PAR attenuation (/mgChl/m^3/m)
inhibcoef= [0.0, 0.0]        # light inhibition
Cquota = [1.8e-11, 1.8e-10]  # Average C quota in cell (mmolC).
Grz_P  = 4000                # phyt.size/Grz_P is the probability to be grazed
dvid_size = 0.5              # relative cell size a cell can start divide
a_dvi = [0.22, 0.22]         # shape parameter for division probability
b_dvi = [2.1, 2.1]           # shape parameter for division probability
death_age = 168.0            # average life time of a cell (7*24 hours)
a_death = 0.05               # shape parameter for natural death
b_death = 0.5                # shape parameter for natural death


Tempref = 293.15   # reference temperature in K
TempAe = -4000.0   # Arrenhius equation
TempCoeff = 0.8    # Arrenhius equation

VNmax = [1.2, 1.2] ./ 86400 # Maximum N uptake rate (mmol N/mmol C/second)
Nqmax = 0.17                # Maximum N quota in cell (mmol N/mmol C)
Nqmin = 0.05                # Minimum N quota in cell (mmol N/mmol C)
KsatN = [0.01, 0.010]       # Half-saturation coeff
VN_b = [0.6, 0.6]        # Shape parameter for size

a_β = 3.1       # scale parameter for metabolic partitioning of biosynthesis
b_β = -3.8      # shape parameter for metabolic pratitioning of biosynthesis
k_mtb = 0.20    # metabolic rate (per hour)

respir_ex= 3.0e-4    # Extra cost of C for biosynthesis
respir_b = 0.13      # Shape parameter for size

grazFracC = 0.7 # Fraction goes into dissolved organic pool
grazFracN = 0.7 # Fraction goes into dissolved organic pool
mortFracC = 0.5 # Fraction goes into dissolved organic pool
mortFracN = 0.5 # Fraction goes into dissolved organic pool

k_sink = 0.01/86400 # Sink rates for agents (m/s)
kDOC   = 1/30/86400 # remineralization rate for DOC, turn over time: a month (per second)
kDON   = 1/30/86400 # remineralization rate for DON, turn over time: a month (per second)
kPOC   = 1/30/86400 # remineralization rate for POC, turn over time: a month (per second)
kPON   = 1/30/86400 # remineralization rate for PON, turn over time: a month (per second)

κh = 0
κv = 0
####################################### end parameters ##################################
