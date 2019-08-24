
## output files 

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

## netcdf

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

## data structures

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

## plotting results

```
julia> using DataFrames, Serialization
julia> output=open(deserialize,"results/output.bin");
julia> plot(output.Population)
```

