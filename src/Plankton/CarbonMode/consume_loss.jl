##### deal with nutrients uptake
function gpu_calc_consume_kernel!(ctsdic, plank, ac, x, y, z, ΔT)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic ctsdic[x[i], y[i], z[i]] += (plank.RS[i] - plank.PS[i]) * ΔT * ac[i]
    end
    return nothing
end
function calc_consume!(ctsdic, plank, ac, x, y, z, ΔT, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_consume_kernel!(ctsdic, plank, ac, x, y, z, ΔT)
    return nothing 
end
function calc_consume!(ctsdic, plank, ac, x, y, z, ΔT, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdic[x[i], y[i], z[i]] += (plank.RS[i] - plank.PS[i]) * ΔT * ac[i]
    end
    return nothing
end

##### deal with grazed or dead individuals
function gpu_calc_loss_kernel!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic ctsdoc[x[i], y[i], z[i]] += plank.Bm[i] * lossFracC * ac[i] * loss[i]
        @inbounds CUDA.@atomic ctspoc[x[i], y[i], z[i]] += plank.Bm[i] * (1.0-lossFracC) * ac[i] * loss[i]
    end
    return nothing
end
function calc_loss!(ctsdoc, ctspoc, plank, ac, x, y, z,
                    loss, lossFracC, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_calc_loss_kernel!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC)
    return nothing 
end
function calc_loss!(ctsdoc, ctspoc, plank, ac, x, y, z, loss, lossFracC, ::CPU)
    for i in 1:size(ac,1)
        @inbounds ctsdoc[x[i], y[i], z[i]] += plank.Bm[i] * lossFracC * ac[i] * loss[i]
        @inbounds ctspoc[x[i], y[i], z[i]] += plank.Bm[i] * (1.0-lossFracC) * ac[i] * loss[i]
    end
    return nothing
end
