module Grids

export AbstractGrid
export RectilinearGrid
export LatLonGrid
export LoadLatLonGrid
export Periodic, Bounded, short_show
export replace_grid_storage
export ΔxC, ΔyC, ΔzC, ΔxF, ΔyF, ΔzF, Ax, Ay, Az, volume 

using Adapt

using PlanktonIndividuals.Architectures

"""
    AbstractGrid{FT, TX, TY, TZ}
Abstract type for grids with elements of type `FT` and topology `{TX, TY, TZ}`.
"""
abstract type AbstractGrid{FT, TX, TY, YZ} end

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