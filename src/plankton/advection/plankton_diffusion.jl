##### calculate diffusivities of each individual
@kernel function plankton_diffusion_kernel!(plank, rnd, κx, κy, κz, ΔT)
    i = @index(Global)
    @inbounds plank.x[i] = plank.x[i] + rnd.x[i] * √(3*ΔT) * 2*√(2*κx) * plank.ac[i]
    @inbounds plank.y[i] = plank.y[i] + rnd.y[i] * √(3*ΔT) * 2*√(2*κy) * plank.ac[i]
    @inbounds plank.z[i] = plank.z[i] + rnd.z[i] * √(3*ΔT) * 2*√(2*κz) * plank.ac[i]
end
function plankton_diffusion!(plank, rnd, κx, κy, κz, ΔT, arch::Architecture)
    kernel! = plankton_diffusion_kernel!(device(arch), 256, (size(plank.ac,1)))
    event = kernel!(plank, rnd, κx, κy, κz, ΔT)
    wait(device(arch), event)

    return nothing
end

plankton_diffusion!(plank, rnd, κ, ΔT, arch) = plankton_diffusion!(plank, rnd, κ, κ, κ, ΔT, arch)

##### generate random numbers for grazing, mortality and division
function gen_rand!(rnd, arch)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    return nothing
end

@kernel function gen_rand_adv_kernel!(rnd)
    i = @index(Global)
    @inbounds rnd.x[i] = rnd.x[i] * 2.0 - 1.0
    @inbounds rnd.y[i] = rnd.y[i] * 2.0 - 1.0
    @inbounds rnd.z[i] = rnd.z[i] * 2.0 - 1.0
end
function gen_rand_adv!(rnd, arch)
    kernel! = gen_rand_adv_kernel!(device(arch), 256, (size(rnd.x,1)))
    event = kernel!(rnd)
    wait(device(arch), event)

    return nothing
end
