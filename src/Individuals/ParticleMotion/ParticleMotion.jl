module ParticleMotion

export particle_advection!, particle_diffusion!
export particle_motion!, colony_motion!
export find_inds!

using KernelAbstractions
using Random

using PlanktonIndividuals.Architectures: device, Architecture, rng_type
using PlanktonIndividuals.Grids
using PlanktonIndividuals: phytoplankton, abiotic_particle

include("interpolation.jl")
include("particle_advection.jl")
include("particle_diffusion.jl")

##### particle advection and diffusion
function particle_motion!(plank::AbstractArray, velos, g::AbstractGrid,
                          vel₀, vel½, vel₁, rnd, κx, κy, κz, ΔT, arch::Architecture)
    particle_advection!(plank, velos, g, vel₀, vel½, vel₁, ΔT, arch)
    particle_diffusion!(plank, rnd, κx, κy, κz, ΔT, g, arch)
    find_inds!(plank, g, arch)
end

##### colony advection and diffusion
function colony_motion!(colony::NamedTuple, velos, g::AbstractGrid, vel₀, vel½, vel₁,
                        rnd, κx, κy, κz, ΔT, arch::Architecture)
    sp1 = colony.sp1
    particle_advection!(sp1.data, velos, g, vel₀, vel½, vel₁, ΔT, arch)
    particle_diffusion!(sp1.data, rnd, κx, κy, κz, ΔT, g, arch)
    find_inds!(sp1.data, g, arch)
    for i in eachindex(colony)[2:end]
        colony[i].data.x  .= sp1.data.x
        colony[i].data.y  .= sp1.data.y
        colony[i].data.z  .= sp1.data.z
        colony[i].data.xi .= sp1.data.xi
        colony[i].data.yi .= sp1.data.yi
        colony[i].data.zi .= sp1.data.zi
    end
end

end