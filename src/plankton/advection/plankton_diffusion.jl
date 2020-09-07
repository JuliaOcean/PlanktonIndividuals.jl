##### calculate diffusivities of each individual
@kernel function calc_diffu_kernel!(phytos, g::Grids, κx, κy, κz, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[1,i] = periodic_domain_x(phytos[1,i] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κx), g)
    @inbounds phytos[2,i] = periodic_domain_y(phytos[2,i] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κy), g)
    @inbounds phytos[3,i] =  bounded_domain_z(phytos[3,i] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κz), g)
end

function plankton_diffusion!(phytos, arch::Architecture, g::Grids, κx, κy, κz, ΔT)
    kernel! = calc_diffu_kernel!(device(arch), 256, (size(phytos,2),))
    event = kernel!(phytos, g, κx, κy, κz, ΔT)
    wait(device(arch), event)
    return nothing
end
plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, ΔT) =
    plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, κ, κ, ΔT)
