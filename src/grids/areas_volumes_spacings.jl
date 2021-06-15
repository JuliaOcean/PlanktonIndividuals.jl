#####
##### Calculate the areas, volumes and spacings of different kind of `AbstractGrid`
##### Currently for RegularRectilinearGrid and RegularLatLonGrid
#####

#####
##### RegularRectilinearGrid
#####
@inline ΔxC(i, j, k, g::RegularRectilinearGrid) = g.Δx
@inline ΔyC(i, j, k, g::RegularRectilinearGrid) = g.Δy
@inline ΔzC(i, j, k, g::RegularRectilinearGrid) = g.Δz

@inline ΔxF(i, j, k, g::RegularRectilinearGrid) = g.Δx
@inline ΔyF(i, j, k, g::RegularRectilinearGrid) = g.Δy
@inline ΔzF(i, j, k, g::RegularRectilinearGrid) = g.Δz

@inline Ax(i, j, k, g::RegularRectilinearGrid) = g.Δy * g.Δz
@inline Ay(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δz
@inline Az(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δy


@inline volume(i, j, k, g::RegularRectilinearGrid) = Az(i, j, k, g) * g.Δz

#####
##### RegularLatLonGrid (degree to meter), halo points included
#####
@inline ΔxC(i, j, k, g::RegularLatLonGrid) = @inbounds g.dxC[i,j]
@inline ΔyC(i, j, k, g::RegularLatLonGrid) = @inbounds g.dyC[i,j]
@inline ΔzC(i, j, k, g::RegularLatLonGrid) = g.Δz

@inline ΔxF(i, j, k, g::RegularLatLonGrid) = @inbounds g.dxF[i,j]
@inline ΔyF(i, j, k, g::RegularLatLonGrid) = @inbounds g.dyF[i,j]
@inline ΔzF(i, j, k, g::RegularLatLonGrid) = ΔzC(i, j, k, g)

@inline Ax(i, j, k, g::RegularLatLonGrid) = @inbounds g.Ax[i,j]
@inline Ay(i, j, k, g::RegularLatLonGrid) = @inbounds g.Ay[i,j]
@inline Az(i, j, k, g::RegularLatLonGrid) = @inbounds g.Az[i,j]

@inline volume(i, j, k, g::RegularLatLonGrid) = @inbounds g.Vol[i,j]

#####
##### VerticallyStretchedLatLonGrid (degree to meter), halo points included
#####
@inline ΔxC(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dxC[i,j]
@inline ΔyC(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dyC[i,j]
@inline ΔzC(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dzC[k]

@inline ΔxF(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dxF[i,j]
@inline ΔyF(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dyF[i,j]
@inline ΔzF(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.dzF[k]

@inline Ax(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.Ax[i,j,k]
@inline Ay(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.Ay[i,j,k]
@inline Az(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.Az[i,j]

@inline volume(i, j, k, g::VerticallyStretchedLatLonGrid) = @inbounds g.Vol[i,j,k]