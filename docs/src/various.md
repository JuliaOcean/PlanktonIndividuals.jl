## Miscellaneous

To build and serve the docs

```
julia make.jl
mkdocs build
mkdocs serve
```

### function inventory

`model_update.jl` is the main program. It runs the time stepping loop + pre- + post-processing.

**Infrastructure functions:**

- `model_includes.jl` include other files (will later be in module file)
- `model_struct.jl` defines `velocity`, `grids`, `nutrient_fields`, and `rem` structs.
- `model_setup.jl` is `setup_agents` and `setup_nutrients`
- `parameters.jl` model parameter values
- `utils.jl` utility functions for IO and operations.

**process functions:**

- `phyt_process.jl` = `daynight `, `PAR_cal `, `PC `, `Nuptake `, `chl_sync `, `divide `, `phyt_update `
- `nutrient_processes.jl` = `compute_nut_biochem`, `compute_source_term`, `nut_update`
- `dst3fl.jl` 3rd order DST Scheme with flux limiting
- `2nd_adv_diffu.jl` right hand side term functions (?)
- `agent_div.jl` = is mostly `agent_move `, `agent_move_1D ` + `double_grid `, `trilinear_itpl `, `simple_itpl `

### model variables

####1) state variables

```
phyts_a
nutrients
```

####2) input variables

- `g` from `grid_offline()` (grid variables)
- `temp, IR` from `read_input()` (before loop)
- `remin` from `rem()` (before loop)
- `velᵇ` from `read_offline_vels()` (in loop)

Notes:

- `IR` and `temp` get cycle through / interpolated inside `phyt_update` (via `IR[trunc(Int,t*ΔT/3600)]`). Not sure about `daynight(t,IR)` in `phyt_update`. 
- `read_offline_vels` receives `trunc(Int,t*ΔT/3600)` as argument.
- `remin` is time-invariant; like `g`.

####3) intermediate variables

- `F` = `compute_nut_biochem(nutrients, remin)`
- `gtr` = `compute_source_term(nutrients, velᵇ, g, F)`
- `nutₜ` = `nut_update(nutrients, consume, g, gtr, ΔT)`

####4) diagnostic variables 

- `dvid_ct`, `graz_ct`, `death_ct` from `phyt_update`
- `gtr`, `nutₜ`, `velᵇ`, `agent_num` from 

### output files 

####1) listing

```
B1.bin		all agents at all time steps for species 1
B2.bin		... species 2
cons_C.txt	total carbon for each time step
cons_DIN.txt	... DIN ...
cons_N.txt	... N ...
grid.bin	grid (input)
IR.bin		irradiance (input)
nutrients/nut.0001.nc (all nitrients for 1 time step)	
nutrients/nut.0002.nc (same...)
output1.bin	average for all agents for each time step for species 1
output2.bin	... 2
output.bin	... for the two species
VD1.bin		vertical profile of the agents opouplation for species 1
VD2.bin		... species 2
```

####2) netcdf output

```
netcdf nut.0001 {
dimensions:
        yC = 1 ;
        xC = 1 ;
        zC = 40 ;
variables:
        double DIC(zC, yC, xC) ;
                DIC:units = "mmolC/m^3" ;
        float yC(yC) ;
                yC:units = "m" ;
                yC:longname = "Locations of the cell centers in the y-direction." ;
        float xC(xC) ;
                xC:units = "m" ;
                xC:longname = "Locations of the cell centers in the x-direction." ;
        float zC(zC) ;
                zC:units = "m" ;
                zC:longname = "Locations of the cell centers in the z-direction." ;
        double DIN(zC, yC, xC) ;
                DIN:units = "mmolN/m^3" ;
        double DOC(zC, yC, xC) ;
                DOC:units = "mmolC/m^3" ;
        double DON(zC, yC, xC) ;
                DON:units = "mmolN/m^3" ;
        double POC(zC, yC, xC) ;
                POC:units = "mmolC/m^3" ;
        double PON(zC, yC, xC) ;
                PON:units = "mmolN/m^3" ;
```

### data structures

```
 output = DataFrame(time=0, 
 gen_ave=mean(B[1].gen), 
 spec_ave = mean(B[1].sp),
 Cq1_ave=mean(B[1].Cq1), 
 Cq2_ave=mean(B[1].Cq2), 
 Nq_ave=mean(B[1].Nq),
 size_ave=mean(B[1].size),
 chl_ave=mean(B[1].chl),
 Population=size(B[1],1),
 dvid=0,
 graz=0,
 death = 0)
```

### plotting results

```
julia> using DataFrames, Serialization
julia> output=open(deserialize,"results/output.bin");
julia> plot(output.Population)
```

