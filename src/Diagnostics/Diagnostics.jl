module Diagnostics

import PlanktonIndividuals: PlanktonDiagnostics

export PlanktonDiagnostics

using KernelAbstractions
using StructArrays

using PlanktonKernels.Architectures: device, Architecture, array_type

using PlanktonIndividuals: Phytoplankton

import Base: show

include("diagnostics_struct.jl")
include("diagnostics_kernels.jl")

end
