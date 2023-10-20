#####
##### record diagnostics of plankton processes at each time step
#####

function gpu_diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:size(ac,1)
        @inbounds CUDA.@atomic diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
    end
end
function diags_proc!(diags_proc, proc, ac, x, y, z, ::GPU)
    @cuda threads=256 blocks=ceil(Int, size(ac,1)/256) gpu_diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    return nothing 
end
function diags_proc!(diags_proc, proc, ac, x, y, z, ::CPU)
    for i in 1:size(ac,1)
        @inbounds diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
    end
end

function diags_spcs!(diags_sp, plank, ac, x, y, z, mode::AbstractMode, arch::Architecture)
    diags = (:PS, :resp, :Bm, :Chl, :Th)
    if isa(mode, CarbonMode)
        diags = (:PS, :RP, :RS, :TD, :Bm, :Bd, :Chl)
    elseif isa(mode, QuotaMode)
        diags = (:PS, :BS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :exu, :Bm, :Cq, :Nq, :Pq, :Chl)
    elseif isa(mode, MacroMolecularMode)
        diags = (:PS, :VDOC, :VHN4, :VNO3, :VPO4, :S_PRO, :S_DNA, :S_RNA, :resp, :ρChl, :CH, :NST, :PST, :PRO, :DNA, :RNA, :Chl)
    end

    for diag in keys(diags_sp)
        if diag in (:num, :graz, :mort, :dvid)
            nothing
        elseif diag in diags
            diags_proc!(diags_sp[diag], getproperty(plank, diag), ac, x, y, z, arch)
        end
    end
end
