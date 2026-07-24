##### grazing and grazing loss
function grazing!(plank, arch::Architecture, plk, p)
    ##### calculate grazing loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data, plk.DFe.data, plk.PFe_bio.data,
               plank, plank.ac, plank.xi, plank.yi, plank.zi, plank.graz, 
               p.grazFracC, p.grazFracN, p.grazFracP, p.grazFracFe, p, arch)
    
    ##### inactivate grazed individuals
    inactivate!(plank, plank.graz)

    return nothing
end

##### mortality and mortality loss
function mortality!(plank, arch::Architecture, plk, p)
    ##### calculate mortality loss
    calc_loss!(plk.DOC.data, plk.POC.data, plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data, plk.DFe.data, plk.PFe_bio.data,
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
        @inbounds plank.PRO_R[idx[i]]  = plank.PRO_R[i]
        @inbounds plank.PRO_Mc[idx[i]] = plank.PRO_Mc[i]
        @inbounds plank.PRO_Mn[idx[i]] = plank.PRO_Mn[i]
        @inbounds plank.PRO_Tn[idx[i]] = plank.PRO_Tn[i]
        @inbounds plank.PRO_Tp[idx[i]] = plank.PRO_Tp[i]
        @inbounds plank.PRO_Tfe[idx[i]] = plank.PRO_Tfe[i]
        @inbounds plank.PRO_RS[idx[i]]  = plank.PRO_RS[i]
        @inbounds plank.PRO_P[idx[i]]  = plank.PRO_P[i]
        @inbounds plank.DNA[idx[i]]  = plank.DNA[i]
        @inbounds plank.RNA[idx[i]]  = plank.RNA[i]
        @inbounds plank.CH[idx[i]]   = plank.CH[i]
        @inbounds plank.qNH3[idx[i]] = plank.qNH3[i]
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
    @inbounds plank.PRO_R[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_Mc[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_Mn[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_Tn[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_Tp[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_Tfe[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.PRO_RS[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.PRO_P[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0
    @inbounds plank.Chl[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.DNA[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.RNA[i] *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.CH[i]  *= (2.0f0 - plank.dvid[i]) / 2.0f0 
    @inbounds plank.qNH3[i]*= (2.0f0 - plank.dvid[i]) / 2.0f0
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
