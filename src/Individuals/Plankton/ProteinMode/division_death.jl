##### grazing and grazing loss
function grazing!(plank, arch::Architecture, plk, p)
    ##### calculate grazing loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, 
               plk.DOP.data, plk.POP.data, plk.DFe.data, plk.PFe_bio.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.grazFracC, p.grazFracN, p.grazFracP, p.grazFracFe, p, arch)
    
    ##### inactivate grazed individuals
    inactivate!(plank, plank.graz)

    return nothing
end

##### mortality and mortality loss
function mortality!(plank, arch::Architecture, plk, p)
    ##### calculate mortality loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, 
               plk.DOP.data, plk.POP.data, plk.DFe.data, plk.PFe_bio.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.mort, 
               p.mortFracC, p.mortFracN, p.mortFracP, p.mortFracFe, p, arch)
    
    ##### inactivate dead individuals
    inactivate!(plank, plank.mort)

    return nothing
end

##### copy ready to divide individuals to inactive rows
@kernel function copy_daughter_individuals_kernel!(plank, con, idx)
    i = @index(Global, Linear)
    if (con[i] == 1.0f0) && (idx[i] ≠ 0)
        @inbounds plank.x[idx[i]]    = plank.x[i]
        @inbounds plank.y[idx[i]]    = plank.y[i]
        @inbounds plank.z[idx[i]]    = plank.z[i]
        @inbounds plank.PRO_RB[idx[i]] = plank.PRO_RB[i]
        @inbounds plank.PRO_MC[idx[i]] = plank.PRO_MC[i]
        @inbounds plank.PRO_MN[idx[i]] = plank.PRO_MN[i]
        @inbounds plank.PRO_TN[idx[i]] = plank.PRO_TN[i]
        @inbounds plank.PRO_TP[idx[i]] = plank.PRO_TP[i]
        @inbounds plank.PRO_TFe[idx[i]]= plank.PRO_TFe[i]
        @inbounds plank.PRO_RS[idx[i]] = plank.PRO_RS[i]
        @inbounds plank.PRO_PS[idx[i]] = plank.PRO_PS[i]
        @inbounds plank.DNA[idx[i]]  = plank.DNA[i]
        @inbounds plank.RNA[idx[i]]  = plank.RNA[i]
        @inbounds plank.CH[idx[i]]   = plank.CH[i]
        @inbounds plank.qNH4[idx[i]] = plank.qNH4[i]
        @inbounds plank.qNO3[idx[i]] = plank.qNO3[i]
        @inbounds plank.qFe[idx[i]]  = plank.qFe[i]
        @inbounds plank.PST[idx[i]]  = plank.PST[i]
        @inbounds plank.Chl[idx[i]]  = plank.Chl[i]
        @inbounds plank.gen[idx[i]]  = plank.gen[i]
        @inbounds plank.ac[idx[i]]   = plank.ac[i]
        @inbounds plank.dvid[idx[i]] = plank.dvid[i]
    end
end
function copy_daughter_individuals!(plank, con, idx::AbstractArray{Int,1}, arch)
    kernel! = copy_daughter_individuals_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, con, idx)
    return nothing
end

##### cell division
@kernel function divide_to_half_kernel!(plank)
    i = @index(Global)
    @inbounds plank.PRO_RB[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_MC[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_MN[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_TN[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_TP[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_TFe[i]*= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_RS[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.PRO_PS[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.Chl[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.DNA[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.RNA[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.CH[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qNH4[i]*= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.qNO3[i]*= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.qFe[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.PST[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.gen[i] += plank.dvid[i]
    @inbounds plank.age[i] *= (1.0f0 - plank.dvid[i])
end
function divide_to_half!(plank, arch)
    kernel! = divide_to_half_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank)
    return nothing
end
function divide!(plank, trs, deactive_ind, arch::Architecture)
    accumulate!(+, trs.idc, plank.dvid)
    trs.idc_int .= unsafe_trunc.(Int, trs.idc)
    get_tind!(plank.idx, plank.dvid, trs.idc_int, deactive_ind, arch)
    copy_daughter_individuals!(plank, plank.dvid, plank.idx, arch)
    divide_to_half!(plank, arch)
    plank.idx .= 0
    return nothing
end
