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

@inline AxC(i, j, k, g::RegularRectilinearGrid) = g.Δy * g.Δz
@inline AyC(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δz
@inline AzC(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δy

@inline AxF(i, j, k, g::RegularRectilinearGrid) = g.Δy * g.Δz
@inline AyF(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δz
@inline AzF(i, j, k, g::RegularRectilinearGrid) = g.Δx * g.Δy

@inline volume(i, j, k, g::RegularRectilinearGrid) = AzF(i, j, k, g) * g.Δz

#####
##### RegularLatLonGrid (degree to meter)
#####
@inline ΔxC(i, j, k, g::RegularLatLonGrid) = @inbounds g.radius * cos(g.yC[j+g.Hy]*π/180) * deg2rad(g.Δx)
@inline ΔyC(i, j, k, g::RegularLatLonGrid) = g.radius * deg2rad(g.Δy)
@inline ΔzC(i, j, k, g::RegularLatLonGrid) = g.Δz

@inline ΔxF(i, j, k, g::RegularLatLonGrid) = ΔxC(i, j, k, g)
@inline ΔyF(i, j, k, g::RegularLatLonGrid) = ΔyC(i, j, k, g)
@inline ΔzF(i, j, k, g::RegularLatLonGrid) = ΔzC(i, j, k, g)

@inline AxC(i, j, k, g::RegularLatLonGrid) = ΔyF(i, j, k, g) * g.Δz
@inline AyC(i, j, k, g::RegularLatLonGrid) = ΔxF(i, j, k, g) * g.Δz
@inline AzC(i, j, k, g::RegularLatLonGrid) = @inbounds g.radius^2 * deg2rad(g.Δx) * (sin(g.yF[j+g.Hy]*π/180) - 
                                                                                     sin(g.yF[j-1+g.Hy]*π/180))

# ΔxF at faces in y direction
@inline ΔxFFA(i, j, k, g::RegularLatLonGrid) = @inbounds g.radius * cos(g.yF[j+g.Hy]*π/180) * deg2rad(g.Δx)

@inline AxF(i, j, k, g::RegularLatLonGrid) = AxF(i, j, k, g)
@inline AyF(i, j, k, g::RegularLatLonGrid) = ΔxFFA(i, j, k, g) *  g.Δz
@inline AzF(i, j, k, g::RegularLatLonGrid) = AzF(i, j, k, g)

@inline volume(i, j, k, g::RegularLatLonGrid) = AzF(i, j, k, g) * g.Δz