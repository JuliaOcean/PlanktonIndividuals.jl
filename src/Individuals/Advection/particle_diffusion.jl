##### calculate diffusivities of each individual
@kernel function calc_diffusion_kernel!(particle, rnd, κx, κy, κz, ΔT, g::AbstractGrid)
    i = @index(Global)
    @inbounds particle.x[i] = particle.x[i] + rnd.x[i] * √(1*κx*ΔT) * particle.ac[i] / 
                                              ΔxC(particle.xi, particle.yi, particle.zi, g)
    @inbounds particle.y[i] = particle.y[i] + rnd.y[i] * √(1*κy*ΔT) * particle.ac[i] / 
                                              ΔyC(particle.xi, particle.yi, particle.zi, g)
    @inbounds particle.z[i] = particle.z[i] + rnd.z[i] * √(1*κz*ΔT) * particle.ac[i] / 
                                              ΔzC(particle.xi, particle.yi, particle.zi, g)
end
function calc_diffusion!(particle, rnd, κx, κy, κz, ΔT, g::AbstractGrid, arch::Architecture)
    kernel! = calc_diffusion_kernel!(device(arch), 256, (size(particle.ac,1)))
    kernel!(particle, rnd, κx, κy, κz, ΔT, g)

    return nothing
end

function particle_diffusion!(particle, rnd, κx, κy, κz, ΔT, g::AbstractGrid, arch::Architecture)
    ##### generate random numbers Normal(0,1)
    randn!(rng_type(arch), rnd.x)
    randn!(rng_type(arch), rnd.y)
    randn!(rng_type(arch), rnd.z)

    ##### calculate diffusion
    calc_diffusion!(particle, rnd, κx, κy, κz, ΔT, arch)

    ##### keep individuals in the domain
    particle_boundaries!(particle, particle.ac, g, arch)
end
