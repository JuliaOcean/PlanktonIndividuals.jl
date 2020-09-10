##### calculate diffusivities of each individual
@kernel function calc_diffu_kernel!(phytos, rand_nums, κx, κy, κz, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[i,1] = phytos[i,1] + rand_nums[i,1] * ΔT * κx
    @inbounds phytos[i,2] = phytos[i,2] + rand_nums[i,2] * ΔT * κy
    @inbounds phytos[i,3] = phytos[i,3] + rand_nums[i,3] * ΔT * κz
end

function calc_diffusion!(phytos, arch::Architecture, κx, κy, κz, ΔT)
    rands = rand(Uniform(-1.0,1.0), size(phytos,1), 3) |> array_type(arch)
    kernel! = calc_diffu_kernel!(device(arch), 256, (size(phytos,1),))
    event = kernel!(phytos, rands, κx, κy, κz, ΔT)
    wait(device(arch), event)
    return nothing
end

function plankton_diffusion!(phytos, arch::Architecture, g::Grids, κx, κy, κz, ΔT)
    calc_diffusion!(phytos, arch, κx, κy, κz, ΔT)
    in_domain!(phytos, arch, g)
end

plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, ΔT) =
    plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, κ, κ, ΔT)
