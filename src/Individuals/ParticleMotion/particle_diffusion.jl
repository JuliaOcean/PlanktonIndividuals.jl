##### calculate diffusivities of each individual
@kernel function calc_diffusion_kernel!(rnd, xi, yi, zi, κx, κy, κz, ΔT, g::AbstractGrid)
    i = @index(Global)
    @inbounds rnd.x[i] = rnd.x[i] * √(κx*ΔT) / ΔxC(xi[i]+g.Hx, yi[i]+g.Hy, zi[i]+g.Hz, g)
    @inbounds rnd.y[i] = rnd.y[i] * √(κy*ΔT) / ΔyC(xi[i]+g.Hx, yi[i]+g.Hy, zi[i]+g.Hz, g)
    @inbounds rnd.z[i] = rnd.z[i] * √(κz*ΔT) / ΔzC(xi[i]+g.Hx, yi[i]+g.Hy, zi[i]+g.Hz, g)
    
end
function calc_diffusion!(rnd, xi, yi, zi, κx, κy, κz, ΔT, g::AbstractGrid, arch::Architecture)
    kernel! = calc_diffusion_kernel!(device(arch), 256, (size(rnd.x,1)))
    kernel!(rnd, xi, yi, zi, κx, κy, κz, ΔT, g)

    return nothing
end

function particle_diffusion!(particle, rnd, κx, κy, κz, ΔT, g::AbstractGrid, arch::Architecture)
    ##### generate random numbers (0,1)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    rnd.x .= (rnd.x .- 0.5f0) .* 2.0f0
    rnd.y .= (rnd.y .- 0.5f0) .* 2.0f0
    rnd.z .= (rnd.z .- 0.5f0) .* 2.0f0
    ##### calculate diffusion
    calc_diffusion!(rnd, particle.xi, particle.yi, particle.zi, κx, κy, κz, ΔT, g, arch)
    particle.x .= particle.x + rnd.x .* particle.ac
    particle.y .= particle.y + rnd.y .* particle.ac
    particle.z .= particle.z + rnd.z .* particle.ac

    ##### keep individuals in the domain
    particle_boundaries!(particle, particle.ac, g, arch)
end
