module Grids

export AbstractGrid
export RegularRectilinearGrid
export RegularLatLonGrid
export VerticallyStretchedLatLonGrid, LoadVerticallyStretchedLatLonGrid
export Periodic, Bounded
export replace_grid_storage
export ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume 

using CUDA
using Adapt

using PlanktonIndividuals.Architectures


abstract type AbstractGrid{TX, TY, YZ} end
abstract type AbstractTopology end
struct Periodic <: AbstractTopology end
struct Bounded <: AbstractTopology end

import Base: show

include("regular_rectilinear_grid.jl")
include("regular_lat_lon_grid.jl")
include("vertically_stretched_lat_lon_grid.jl")
include("utils.jl")
include("areas_volumes_spacings.jl")

end