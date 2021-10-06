module Carbon

using KernelAbstractions: @kernel, @index, Event, MultiEvent, wait
using KernelAbstractions.Extras.LoopInfo: @unroll
using CUDA
using StructArrays
using Random
using LinearAlgebra

using PlanktonIndividuals.Architectures: device, Architecture, GPU, CPU, rng_type, array_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals.Output

using PlanktonIndividuals: AbstractMode, CarbonMode, QuotaMode, MacroMolecularMode

end