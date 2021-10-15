# Model Simulation

A `PlanktonSimulation` includes a `PlanktonModel` and its time step `ΔT` and number of time steps `nΔT`.
It will time step the `PlanktonModel` by calling [`update!`](@ref).

```@docs
PlanktonSimulation
```

## Model Diagnostics

```@docs
PlanktonDiagnostics
```

Model diagnostics are specified by `tracer` (for tracers) and `plankton` (for individuals).
Diagnostics for individuals are aggregated into gridded fields.

A full list of available diagnostics is provided below:

```julia
tracer = (:PAR, # photosynthetically active radiation
          :DIC, # dissolved inorganic carbon
          :NH4, # ammonia
          :NO3, # nitrate
          :PO4, # phosphate
          :DOC, # dissolved organic carbon
          :DON, # dissolved organic nitrogen
          :DOP, # dissolved organic phosphorus
          :POC, # particulate organic carbon
          :PON, # particulate organic nitrogen
          :POP  # particulate organic phosphorus
         )

plankton = (:num,  # number of individuals
            :graz, # number of grazed individuals
            :mort, # number of died individuals
            :dvid, # number of divided individuals
            :PS,   # photosynthesis rate
            :BS,   # biosynthesis rate
            :VDOC, # DOC uptake rate
            :VHN4, # NH4 uptake rate
            :VNO3, # NO3 uptake rate
            :VPO4, # PO4 uptake rate
            :resp, # respiration rate
            :exu,  # exudation rate
            :Bm,   # functional biomass
            :Cq,   # Carbon pool
            :Nq,   # Nitrogen pool
            :Pq,   # Phosphorus pool
            :Chl   # Chla
           )
```

## Output Writer

```@docs
PlanktonOutputWriter
```

The model currently has two types of outputs which are both saved in `JLD2` files.

First, the state of all `individuals`
at each time step of a `PlanktonSimulation` gets saved in a file named `plankton_prefix*".jld2"`, default: `plankton.jld2`.

Second, for diagnostics, `individuals` at each time step will be aggregated into tracer fields.
The keyword argument `time_interval` in `PlanktonDiagnostics` specifies the time window that the diagnostics are averaged.
Only diagnostics specified by `tracer` and `plankton` in `PlanktonDiagnostics` will be saved.
All the diagnostics of a `PlanktonSimulation` will be saved in a file named `diags_prefix*".jld2"`, default: `diags.jld2`.
