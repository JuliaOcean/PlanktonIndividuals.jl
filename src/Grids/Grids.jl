module Grids

export AbstractGrid
export RectilinearGrid
export LatLonGrid
export LoadLatLonGrid
export Periodic, Bounded, short_show
export replace_grid_storage
export ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume 

using CUDA
using Adapt

using PlanktonIndividuals.Architectures

"""
    AbstractGrid{TX, TY, TZ}
Abstract type for grids with elements of type `Float64` and topology `{TX, TY, TZ}`.
"""
abstract type AbstractGrid{TX, TY, YZ} end

"""
    AbstractTopology
Abstract type for grid topologies.
"""
abstract type AbstractTopology end

"""
    Periodic
Grid topology for periodic dimensions.
"""
struct Periodic <: AbstractTopology end

"""
    Bounded
Grid topology for bounded dimensions.
"""
struct Bounded <: AbstractTopology end

import Base: show

include("rectilinear_grid.jl")
include("lat_lon_grid.jl")
include("utils.jl")
include("areas_volumes_spacings.jl")

end