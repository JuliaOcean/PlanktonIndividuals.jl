##### calculate diffusivities of each individual
@kernel function calc_diffu_kernel!(plank, rand_nums, κx, κy, κz, ΔT)
    i = @index(Global, Linear)
    @inbounds plank[i,1] = plank[i,1] + rand_nums[i,1] * √(3*ΔT) * 2*√(2*κx)
    @inbounds plank[i,2] = plank[i,2] + rand_nums[i,2] * √(3*ΔT) * 2*√(2*κx)
    @inbounds plank[i,3] = plank[i,3] + rand_nums[i,3] * √(3*ΔT) * 2*√(2*κx)
end

function calc_diffusion!(plank, arch::Architecture, κx, κy, κz, ΔT)
    rands = rand(Uniform(-1.0,1.0), size(plank,1), 3) |> array_type(arch)
    kernel! = calc_diffu_kernel!(device(arch), 256, (size(plank,1),))
    event = kernel!(plank, rands, κx, κy, κz, ΔT)
    wait(device(arch), event)
    return nothing
end

function plankton_diffusion!(plank, arch::Architecture, g::Grids, κx, κy, κz, ΔT)
    calc_diffusion!(plank, arch, κx, κy, κz, ΔT)
    in_domain!(plank, arch, g)
end

plankton_diffusion!(plank, arch::Architecture, g::Grids, κ, ΔT) =
    plankton_diffusion!(plank, arch::Architecture, g::Grids, κ, κ, κ, ΔT)
