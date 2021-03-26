##### calculate diffusivities of each individual
@kernel function calc_diffusion_kernel!(plank, rnd, κx, κy, κz, ΔT)
    i = @index(Global)
    @inbounds plank.x[i] = plank.x[i] + rnd.x[i] * √(1*κx*ΔT) * plank.ac[i]
    @inbounds plank.y[i] = plank.y[i] + rnd.y[i] * √(1*κy*ΔT) * plank.ac[i]
    # @inbounds plank.z[i] = plank.z[i] + rnd.z[i] * √(1*κz*ΔT) * plank.ac[i]
end
function calc_diffusion!(plank, rnd, κx, κy, κz, ΔT, arch::Architecture)
    kernel! = calc_diffusion_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, rnd, κx, κy, κz, ΔT)
    wait(device(arch), event)

    return nothing
end

function plankton_diffusion!(plank, rnd, κx, κy, κz, ΔT, g::RegularRectilinearGrid, arch::Architecture)
    ##### generate random numbers for grazing, mortality and division (0,1)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    ##### shift the random numbers to (-1,1)
    @inbounds rnd.x .= rnd.x .* 2.0 .- 1.0
    @inbounds rnd.y .= rnd.y .* 2.0 .- 1.0
    @inbounds rnd.z .= rnd.z .* 2.0 .- 1.0

    ##### calculate diffusion
    calc_diffusion!(plank, rnd, κx, κy, κz, ΔT, arch)

    ##### keep individuals in the domain
    periodic_domain!(plank, plank.ac, g, arch)
end

plankton_diffusion!(plank, rnd, κ, ΔT, g, arch) = plankton_diffusion!(plank, rnd, κ, κ, κ, ΔT, g, arch)