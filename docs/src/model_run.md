# Model Simulation

A `PI_simulation` includes a `PI_model` and its time steps `nΔT`, `ΔT`.
It will time step the `PI_model` by calling [`update!'](@ref).

```@docs
PI_simulation
```

## Diagnostics

Model diagnostics are specified by `diag_ntrs` (for tracers) and `diag_nprocs` (for individuals).
Diagnostics for individuals are aggregated into tracer fields.

A full list of available diagnostics are provided below:

```julia
diag_ntrs = (:PAR, :DIC, :NH4, :NO3, :PO4, :DOC, :DON, :DOP, :POC, :PON, :POP)

diag_nprocs = (:num,  # number of individuals
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
               :chl   # Chla
               )
```

## Output

Currently, we only support two types of output which are both save in `JLD2` files.
The first type of output is for `individual`s. The current state of all the `individuals`
at each time step of a `PI_simulation` will be saved in a single file named `individuals.jld2`.
An example structure of `individuals.jld2` is shown below.

```julia
julia> jldopen("results/individuals.jld2") # only the location and cell size is saved for now
JLDFile /home/zhenwu/PI_GPU/results/individuals.jld2 (read-only)
 ├─📂 0000000060
 │  └─📂 sp1
 │     ├─🔢 x
 │     ├─🔢 y
 │     ├─🔢 z
 │     └─🔢 Sz
 └─📂 0000000120
    └─📂 sp1
       ├─🔢 x
       ├─🔢 y
       ├─🔢 z
       └─🔢 Sz
```

The second type of output is for diagnostics. `individual`s at each time step will be aggregated into tracer fields.
The frequency of diagnostics is specified by `diag_freq` in `PI_simulation`.
Only diagnostics specified by `diag_ntrs` and `diag_nprocs` in `PI_Model` will be saved.
All the diagnostics of a `PI_simulation` will be saved in a single file named `diags.jld2`.
An example structure of `diags.jld2` is shown below.

```julia
julia> jldopen("results/diags.jld2")
JLDFile /home/zhenwu/PI_GPU/results/diags.jld2 (read-only)
 ├─📂 0000000060
 │  ├─📂 nut
 │  │  ├─🔢 PAR
 │  │  ├─🔢 DOC
 │  │  ├─🔢 NH4
 │  │  └─🔢 NO3
 │  └─📂 sp1
 │     ├─🔢 num
 │     ├─🔢 graz
 │     ├─🔢 mort
 │     └─🔢 dvid
 └─📂 0000000120
    ├─📂 nut
    │  ├─🔢 PAR
    │  ├─🔢 DOC
    │  ├─🔢 NH4
    │  └─🔢 NO3
    └─📂 sp1
       ├─🔢 num
       ├─🔢 graz
       ├─🔢 mort
       └─🔢 dvid
```
