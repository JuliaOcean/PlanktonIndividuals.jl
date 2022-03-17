##### deal with nutrients uptake
function gpu_calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic ctsdic[x[i], y[i], z[i]] += (plank.resp[i] - plank.PS[i]) * ΔT * ac[i]
        @inbounds CUDA.@atomic ctsdoc[x[i], y[i], z[i]] += (plank.exu[i] - plank.VDOC[i]) * ΔT * ac[i]
        @inbounds CUDA.@atomic ctsnh4[x[i], y[i], z[i]] += -plank.VNH4[i] * ΔT * ac[i]
        @inbounds CUDA.@atomic ctsno3[x[i], y[i], z[i]] += -plank.VNO3[i] * ΔT * ac[i]
        @inbounds CUDA.@atomic ctspo4[x[i], y[i], z[i]] += -plank.VPO4[i] * ΔT * ac[i]
    end
    return nothing
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_consume_kernel!(ctsdic, ctsdoc, ctsnh4, 
                                                        ctsno3, ctspo4, plank, ac, x, y, z, ΔT)
    return nothing 
end
function calc_consume!(ctsdic, ctsdoc, ctsnh4, ctsno3, ctspo4, plank, ac, x, y, z, ΔT, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdic[x[i], y[i], z[i]] += (plank.resp[i] - plank.PS[i]) * ΔT * ac[i]
        @inbounds ctsdoc[x[i], y[i], z[i]] += (plank.exu[i] - plank.VDOC[i]) * ΔT * ac[i]
        @inbounds ctsnh4[x[i], y[i], z[i]] += -plank.VNH4[i] * ΔT * ac[i]
        @inbounds ctsno3[x[i], y[i], z[i]] += -plank.VNO3[i] * ΔT * ac[i]
        @inbounds ctspo4[x[i], y[i], z[i]] += -plank.VPO4[i] * ΔT * ac[i]
    end
    return nothing
end

##### deal with grazed or dead individuals
function gpu_calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                                   loss, lossFracC, lossFracN, lossFracP, p)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic ctsdoc[x[i], y[i], z[i]] += total_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.CH[i]) * lossFracC * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctsdon[x[i], y[i], z[i]] += total_N_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.NST[i], p) * lossFracN * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctsdop[x[i], y[i], z[i]] += total_P_biomass(plank.DNA[i], plank.RNA[i], plank.PST[i], p) * lossFracP * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctspoc[x[i], y[i], z[i]] += total_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.CH[i]) * (1.0-lossFracC) * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctspon[x[i], y[i], z[i]] += total_N_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.NST[i], p) * (1.0-lossFracN) * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctspop[x[i], y[i], z[i]] += total_P_biomass(plank.DNA[i], plank.RNA[i], plank.PST[i], p) * (1.0-lossFracP) * ac[i] * loss[i]
    end
    return nothing
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                    loss, lossFracC, lossFracN, lossFracP, p, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_loss_kernel!(ctsdoc, ctspoc, ctsdon, ctspon, 
                                                        ctsdop, ctspop, plank, ac, x, y, z, loss, 
                                                        lossFracC, lossFracN, lossFracP, p)
    return nothing 
end
function calc_loss!(ctsdoc, ctspoc, ctsdon, ctspon, ctsdop, ctspop, plank, ac, x, y, z,
                           loss, lossFracC, lossFracN, lossFracP, p, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdoc[x[i], y[i], z[i]] += total_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.CH[i]) * lossFracC * ac[i] * loss[i]
        @inbounds ctsdon[x[i], y[i], z[i]] += total_N_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.NST[i], p) * lossFracN * ac[i] * loss[i]
        @inbounds ctsdop[x[i], y[i], z[i]] += total_P_biomass(plank.DNA[i], plank.RNA[i], plank.PST[i], p) * lossFracP * ac[i] * loss[i]
        @inbounds ctspoc[x[i], y[i], z[i]] += total_C_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.CH[i]) * (1.0-lossFracC) * ac[i] * loss[i]
        @inbounds ctspon[x[i], y[i], z[i]] += total_N_biomass(plank.PRO[i], plank.DNA[i], plank.RNA[i], plank.NST[i], p) * (1.0-lossFracN) * ac[i] * loss[i]
        @inbounds ctspop[x[i], y[i], z[i]] += total_P_biomass(plank.DNA[i], plank.RNA[i], plank.PST[i], p) * (1.0-lossFracP) * ac[i] * loss[i]
    end
    return nothing
end