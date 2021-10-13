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

function diags_spcs!(diags_sp, proc, plank, ac, x, y, z, mode::AbstractMode, arch::Architecture)
    proc_diags = (:PS, :resp)
    plank_diags = (:Bm, :Chl)
    if isa(mode, CarbonMode)
        nothing
    elseif isa(mode, QuotaMode)
        proc_diags = (:PS, :BS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :exu)
        plank_diags = (:Bm, :Cq, :Nq, :Pq, :Chl)
    elseif isa(mode, MacroMolecularMode)
        nothing
    end

    for diag in keys(diags_sp)
        if diag in (:num, :graz, :mort, :dvid)
            nothing
        elseif diag in proc_diags
            diags_proc!(diags_sp[diag], getproperty(proc, diag), ac, x, y, z, arch)
        elseif diag in plank_diags
            diags_proc!(diags_sp[diag], getproperty(plank, diag), ac, x, y, z, arch)
        end
    end
end
