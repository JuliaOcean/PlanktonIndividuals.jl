##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ctsdic, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, R_NC, R_PC)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdic[x[i], y[i], z[i]] -= plank.BS[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsnh4[x[i], y[i], z[i]] -= plank.BS[i] * R_NC * ΔT * ac[i] *
                                                                        (plank.fNH4[i]/max(1.0f-30, plank.fNH4[i] + plank.fNO3[i]))
    @inbounds KernelAbstractions.@atomic ctsno3[x[i], y[i], z[i]] -= plank.BS[i] * R_NC * ΔT * ac[i] *
                                                                        (plank.fNO3[i]/max(1.0f-30, plank.fNH4[i] + plank.fNO3[i]))
    @inbounds KernelAbstractions.@atomic ctspo4[x[i], y[i], z[i]] -= plank.BS[i] * R_PC * ΔT * ac[i]
end
function calc_consume!(ctsdic, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, R_NC, R_PC, arch)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdic, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, R_NC, R_PC)
    return nothing 
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                           loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdoc[x[i], y[i], z[i]] += plank.Bm[i] * lossFracC * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdon[x[i], y[i], z[i]] += plank.Bm[i] * R_NC * lossFracN * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdop[x[i], y[i], z[i]] += plank.Bm[i] * R_PC * lossFracP * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspoc[x[i], y[i], z[i]] += plank.Bm[i] * (1.0f0-lossFracC) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspon[x[i], y[i], z[i]] += plank.Bm[i] * R_NC * (1.0f0-lossFracN) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspop[x[i], y[i], z[i]] += plank.Bm[i] *R_PC * (1.0f0-lossFracP) * ac[i] * loss[i]
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
    loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC, arch)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z, 
            loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
return nothing 
end