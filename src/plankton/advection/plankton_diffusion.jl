##### calculate diffusivities of each individual
# @kernel function calc_diffu_kernel!(plank, rnd, κx, κy, κz, ΔT)
#     i = @index(Global, Linear)
#     if plank[i,58] == 1.0
#         @inbounds plank[i,1] = plank[i,1] + rnd[i,1] * √(3*ΔT) * 2*√(2*κx)
#         @inbounds plank[i,2] = plank[i,2] + rnd[i,2] * √(3*ΔT) * 2*√(2*κx)
#         @inbounds plank[i,3] = plank[i,3] + rnd[i,3] * √(3*ΔT) * 2*√(2*κx)
#     end
# end

# function plankton_diffusion!(plank, rnd, arch::Architecture, κx, κy, κz, ΔT)
#     kernel! = calc_diffu_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank, rnd, κx, κy, κz, ΔT)
#     wait(device(arch), event)
#     return nothing
# end

# plankton_diffusion!(plank, rnd, arch, κ, ΔT) = plankton_diffusion!(plank, rnd, arch, κ, κ, κ, ΔT)

function plankton_diffusion!(plank, rnd, κx, κy, κz, ΔT)
    plank.x .= plank.x .+ rnd.x .* √(3*ΔT) .* (2*√(2*κx)) .* plank.ac
    plank.y .= plank.y .+ rnd.y .* √(3*ΔT) .* (2*√(2*κy)) .* plank.ac
    plank.z .= plank.z .+ rnd.z .* √(3*ΔT) .* (2*√(2*κz)) .* plank.ac

    return nothing
end

plankton_diffusion!(plank, rnd, κ, ΔT) = plankton_diffusion!(plank, rnd, κ, κ, κ, ΔT)

function gen_rand_adv!(rnd, arch)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)
    rnd.x .= rnd.x .* 2.0 .- 1.0
    rnd.y .= rnd.y .* 2.0 .- 1.0
    rnd.z .= rnd.z .* 2.0 .- 1.0
end
