##### calculate diffusivities of each individual
@kernel function calc_diffu_kernel!(phytos, g::Grids, κx, κy, κz, ΔT)
    i = @index(Global, Linear)
    @inbounds phytos[i,1] = periodic_domain_x(phytos[i,1] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κx), g)
    @inbounds phytos[i,2] = periodic_domain_y(phytos[i,2] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κy), g)
    @inbounds phytos[i,3] =  bounded_domain_z(phytos[i,3] + rand(-1.0:0.001:1.0) * √(3*ΔT) * 2*√(2*κz), g)
end

function plankton_diffusion!(phytos, arch::Architecture, g::Grids, κx, κy, κz, ΔT)
    kernel! = calc_diffu_kernel!(device(arch), 256, (size(phytos,1),))
    event = kernel!(phytos, g, κx, κy, κz, ΔT)
    wait(device(arch), event)
    return nothing
end
plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, ΔT) =
    plankton_diffusion!(phytos, arch::Architecture, g::Grids, κ, κ, κ, ΔT)
