#####
##### Calculate the areas, volumes and spacings of different kind of `AbstractGrid`
##### Currently for RectilinearGrid and LatLonGrid
#####

#####
##### RectilinearGrid, halo points included
#####
@inline ΔxC(i, j, k, g::RectilinearGrid) = g.Δx
@inline ΔyC(i, j, k, g::RectilinearGrid) = g.Δy
@inline ΔzC(i, j, k, g::RectilinearGrid) = @inbounds g.dzC[k]

@inline ΔxF(i, j, k, g::RectilinearGrid) = g.Δx
@inline ΔyF(i, j, k, g::RectilinearGrid) = g.Δy
@inline ΔzF(i, j, k, g::RectilinearGrid) = @inbounds g.dzF[k]

@inline Ax(i, j, k, g::RectilinearGrid) = @inbounds g.Δy * g.dzF[k]
@inline Ay(i, j, k, g::RectilinearGrid) = @inbounds g.Δx * g.dzF[k]
@inline Az(i, j, k, g::RectilinearGrid) = g.Δx * g.Δy


@inline volume(i, j, k, g::RectilinearGrid) = @inbounds g.Δx * g.Δy * g.dzF[k]

#####
##### LatLonGrid (degree to meter), halo points included
#####
@inline ΔxC(i, j, k, g::LatLonGrid) = @inbounds g.dxC[i,j]
@inline ΔyC(i, j, k, g::LatLonGrid) = @inbounds g.dyC[i,j]
@inline ΔzC(i, j, k, g::LatLonGrid) = @inbounds g.dzC[k]

@inline ΔxF(i, j, k, g::LatLonGrid) = @inbounds g.dxF[i,j]
@inline ΔyF(i, j, k, g::LatLonGrid) = @inbounds g.dyF[i,j]
@inline ΔzF(i, j, k, g::LatLonGrid) = @inbounds g.dzF[k]

@inline Ax(i, j, k, g::LatLonGrid) = @inbounds g.Ax[i,j,k]
@inline Ay(i, j, k, g::LatLonGrid) = @inbounds g.Ay[i,j,k]
@inline Az(i, j, k, g::LatLonGrid) = @inbounds g.Az[i,j]

@inline volume(i, j, k, g::LatLonGrid) = @inbounds g.Vol[i,j,k]