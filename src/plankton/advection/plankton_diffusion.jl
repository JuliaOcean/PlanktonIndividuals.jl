##### calculate diffusivities of each individual
@kernel function calc_diffu_kernel!(plank, κx, κy, κz, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,1] = plank[i,1] + plank[i,58] * √(3*ΔT) * 2*√(2*κx)
    @inbounds plank[i,2] = plank[i,2] + plank[i,59] * √(3*ΔT) * 2*√(2*κx)
    @inbounds plank[i,3] = plank[i,3] + plank[i,60] * √(3*ΔT) * 2*√(2*κx)
end

function calc_diffusion!(plank, arch::Architecture, κx, κy, κz, ΔT)
    kernel! = calc_diffu_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, κx, κy, κz, ΔT)
    wait(device(arch), event)
    return nothing
end

# @kernel function gen_rand_kernel!(plank)
#     i = @index(Global, Linear)
#     @inbounds plank[i,58] = rand(Uniform(-1.0,1.0))
#     @inbounds plank[i,59] = rand(Uniform(-1.0,1.0))
#     @inbounds plank[i,60] = rand(Uniform(-1.0,1.0))
# end
# function gen_rand!(plank, arch::Architecture)
#     kernel! = gen_rand_kernel!(device(arch), 256, (size(plank,1),))
#     event = kernel!(plank)
#     wait(device(arch), event)
#     return nothing
# end

function plankton_diffusion!(plank, arch::Architecture, g::Grids, κx, κy, κz, ΔT)
    # gen_rand!(plank, arch)
    plank[:,58:60] .= rand(Uniform(-1.0,1.0), size(plank,1), 3) |> array_type(arch)
    calc_diffusion!(plank, arch, κx, κy, κz, ΔT)
    in_domain!(plank, arch, g)
end

plankton_diffusion!(plank, arch::Architecture, g::Grids, κ, ΔT) =
    plankton_diffusion!(plank, arch::Architecture, g::Grids, κ, κ, κ, ΔT)
