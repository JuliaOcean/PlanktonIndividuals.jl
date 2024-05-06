##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ctsdic, plank, ac, x, y, z, ΔT)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdic[x[i], y[i], z[i]] += (plank.RS[i] - plank.PS[i]) * ΔT * ac[i]
end
function calc_consume!(ctsdic, plank, ac, x, y, z, ΔT, arch)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdic, plank, ac, x, y, z, ΔT)
    return nothing 
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdoc[x[i], y[i], z[i]] += plank.Bm[i] * lossFracC * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspoc[x[i], y[i], z[i]] += plank.Bm[i] * (1.0-lossFracC) * ac[i] * loss[i]
end
function calc_loss!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC, arch)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC)
return nothing 
end