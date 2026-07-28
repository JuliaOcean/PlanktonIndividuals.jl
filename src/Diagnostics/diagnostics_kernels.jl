#####
##### record diagnostics of particle processes at each time step
#####
@kernel function diags_proc_kernel!(diags_proc, proc, ac, x, y, z)
    i = @index(Global)
    @inbounds KernelAbstractions.@atomic diags_proc[x[i], y[i], z[i]] += proc[i] * ac[i]
end
function diags_proc!(diags_proc, proc, ac, x, y, z, arch)
    kernel! = diags_proc_kernel!(device(arch), 256, (size(ac,1)))
    kernel!(diags_proc, proc, ac, x, y, z)
    return nothing 
end

function diags_spcs!(diags_sp, plank::phytoplankton, ac, x, y, z, mode::AbstractMode, arch::Architecture)
    diags = (:PS, :resp, :Bm, :Chl)
    if isa(mode, CarbonMode)
        diags = (:PS, :BS, :RP, :RS, :TD, :Bm, :Bd, :Chl, :ptc)
    elseif isa(mode, QuotaMode)
        diags = (:PS, :BS, :VDOC, :VNH4, :VNO3, :VPO4, :resp, :exu, :Bm, :Cq, :Nq, :Pq, :Chl, :ptc)
    elseif isa(mode, MacroMolecularMode)
        diags = (:PS, :VDOC, :VHN4, :VNO3, :VPO4, :S_PRO, :S_DNA, :S_RNA, :resp, :ρChl, 
                 :CH, :NST, :PST, :PRO, :DNA, :RNA, :Chl, :ptc)
    elseif isa(mode, IronEnergyMode)
        diags = (:PS, :CF, :ECF, :RS, :ERS, :NR, :ENR, :NF, :ENF, :BS, :VNH4, :VNO3, :VPO4, :VFe, 
                 :PS2ST, :ST2PS, :NR2ST, :ST2NR, :NF2ST, :ST2NF, :Bm, :exEn, :CH, :qNO3, :qNH4, 
                 :qP, :qFe, :qFePS, :qFeNR, :qFeNF, :Chl, :tdark, :ptc)
    elseif isa(mode, ProteinMode)
        diags = (:PS, :CF, :ECF, :RS, :ERS, :NR, :ENR, :NF, :ENF,
                 :VNH4, :VNO3, :VPO4, :VFe, :VDOC,
                 :PRO_RB, :PRO_MC, :PRO_MN, :PRO_TN, :PRO_TP, :PRO_TFe, :PRO_RS, :PRO_PS,
                 :DNA, :RNA, :CH, :Chl, :qNH4, :qNO3, :qFe, :PST,
                 :SP_RB, :SP_MC, :SP_MN, :SP_TN, :SP_TP, :SP_TFe, :SP_RS, :SP_PS, :SDNA, :SRNA,
                 :ESP_RB, :ESP_MC, :ESP_MN, :ESP_TN, :ESP_TP, :ESP_TFe, :ESP_RS, :ESP_PS, :EDNA, :ERNA,
                 :exE_RS, :exE_PS,
                 :DP_RB, :DP_PS, :DP_MC, :DChl,
                 :ρChl, :exu, :ptc)
    end

    for diag in keys(diags_sp)
        if diag in (:num, :graz, :mort, :dvid, :ptc)
            nothing
        elseif diag in diags
            diags_proc!(diags_sp[diag], getproperty(plank.data, diag), ac, x, y, z, arch)
        end
    end
end