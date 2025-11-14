##### deal with nutrients uptake
@kernel function calc_consume_kernel!(ctsdic, ctsnh4, ctsno3, ctspo4, ctsdfe, csto2,
                                        plank, ac, x, y, z, ΔT)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdic[x[i], y[i], z[i]] += (plank.RS[i] - plank.PS[i]) * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsnh4[x[i], y[i], z[i]] += -plank.VNH4[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsno3[x[i], y[i], z[i]] += -plank.VNO3[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctspo4[x[i], y[i], z[i]] += -plank.VPO4[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsdfe[x[i], y[i], z[i]] += -plank.VFe[i]  * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctso2[x[i], y[i], z[i]]  += -plank.VO2[i]  * ΔT * ac[i]
end
function calc_consume!(ctsdic, ctsnh4, ctsno3, ctspo4, ctsdfe, csto2, plank, ac, x, y, z, ΔT, arch)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdic, ctsnh4, ctsno3, ctspo4, ctsdfe, csto2, plank, ac, x, y, z, ΔT)
    return nothing 
end

##### deal with grazed or dead individuals
@kernel function calc_loss_kernel!(ctsnh4, ctsno3, ctsdoc, ctspoc, ctsdon, ctspon, 
                                    ctsdop, ctspop, ctsdfe, ctspfe_bio, csto2,
                                    plank, ac, x, y, z, loss, 
                                    lossFracC, lossFracN, lossFracP, lossFracFe, R_NC, R_PC)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsnh4[x[i], y[i], z[i]] += plank.qNH4[i] * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsno3[x[i], y[i], z[i]] += plank.qNO3[i] * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.CH[i]) * lossFracC * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdon[x[i], y[i], z[i]] += plank.Bm[i]*R_NC * lossFracN * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.qP[i]) * lossFracP * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdfe[x[i], y[i], z[i]] += (plank.qFe[i] + plank.qFePS[i] + plank.qFeNF[i] + plank.qFeNR[i]) * lossFracFe * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.CH[i]) * (1.0f0-lossFracC) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspon[x[i], y[i], z[i]] += plank.Bm[i]*R_NC * (1.0f0-lossFracN) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.qP[i]) * (1.0f0-lossFracP) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspfe_bio[x[i], y[i], z[i]] += (plank.qFe[i] + plank.qFePS[i] + plank.qFeNF[i] + plank.qFeNR[i]) * (1.0f0-lossFracFe) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctso2[x[i], y[i], z[i]]  += plank.qO2[i] * ac[i] * loss[i]
end
function calc_loss!(ctsnh4, ctsno3, ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, ctsdfe, ctspfe_bio, csto2,
    plank, ac, x, y, z, loss, lossFracC, lossFracN, lossFracP, lossFracFe, R_NC, R_PC, arch)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsnh4, ctsno3, ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, ctsdfe, ctspfe_bio, csto2,
            plank, ac, x, y, z, loss, lossFracC, lossFracN, lossFracP, lossFracFe, R_NC, R_PC)
return nothing 
end