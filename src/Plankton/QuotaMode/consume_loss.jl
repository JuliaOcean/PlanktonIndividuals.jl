##### deal with nutrients uptake
function calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdic[x[i], y[i], z[i]] += (plank.resp[i] - plank.PS[i]) * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsdoc[x[i], y[i], z[i]] += (plank.exu[i] - plank.VDOC[i]) * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsnh4[x[i], y[i], z[i]] += -plank.VNH4[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctsno3[x[i], y[i], z[i]] += -plank.VNO3[i] * ΔT * ac[i]
    @inbounds KernelAbstractions.@atomic ctspo4[x[i], y[i], z[i]] += -plank.VPO4[i] * ΔT * ac[i]
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, arch)
    kernel! = calc_consume_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT)
    return nothing 
end

##### deal with grazed or dead individuals
function calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                           loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic ctsdoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * lossFracC * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * lossFracN * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctsdop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * lossFracP * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspoc[x[i], y[i], z[i]] += (plank.Bm[i] + plank.Cq[i]) * (1.0-lossFracC) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspon[x[i], y[i], z[i]] += (plank.Bm[i]*R_NC + plank.Nq[i]) * (1.0-lossFracN) * ac[i] * loss[i]
    @inbounds KernelAbstractions.@atomic ctspop[x[i], y[i], z[i]] += (plank.Bm[i]*R_PC + plank.Pq[i]) * (1.0-lossFracP) * ac[i] * loss[i]
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
    loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC, arch)
    kernel! = calc_loss_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z, 
            loss, lossFracC, lossFracN, lossFracP, R_NC, R_PC)
return nothing 
end