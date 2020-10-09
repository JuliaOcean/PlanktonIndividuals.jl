##### copy active individuals to model.timestepper.tmp
@kernel function copyto_tmp_kernel!(plank, tmp, con, idx, b::Bool)
    i,j = @index(Global, NTuple)
    if con[i] == 1.0
        @inbounds tmp[idx[i],j] = copy(plank[i,j])
        @inbounds plank[i,j] = plank[i,j] * b # (de)activate individual
    end
end
function copyto_tmp!(plank, tmp, con, idx::AbstractArray{Int64,1}, b::Bool, arch::Architecture)
    kernel! = copyto_tmp_kernel!(device(arch), (16,16), (size(plank,1),60))
    event = kernel!(plank, tmp, con, idx, b)
    wait(device(arch), event)
    return nothing
end

##### grazing and grazing loss
function grazing!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank[:,59] .= 0.0
    plank[:,59] .= cumsum(plank[:,31])
    ##### copy grazed individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank[:,31], Int.(plank[:,59]), false, arch)
    ##### calculate grazing loss
    calc_loss!(tmp, Int.(tmp[:,13:15]), arch, plk.DOC.data, plk.POC.data,
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.grazFracC, p.grazFracN, p.grazFracP, p.R_NC, p.R_PC)
end

##### mortality and mortality loss
function mortality!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank[:,59] .= 0.0
    plank[:,59] .= cumsum(plank[:,32])
    ##### copy grazed individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank[:,32], Int.(plank[:,59]), false, arch)
    ##### calculate grazing loss
    calc_loss!(tmp, Int.(tmp[:,13:15]), arch, plk.DOC.data, plk.POC.data,
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.mortFracC, p.mortFracN, p.mortFracP, p.R_NC, p.R_PC)
end

##### cell division
function divide_copy!(tmp, arch::Architecture, tmp_num::Int64)
    ##### calculate index for timestepper.tmp
    tmp[:,59] .= 0.0
    tmp[:,59] .= cumsum(tmp[:,33])
    tmp[:,59] .= tmp[:,59] .+ tmp_num
    copyto_tmp!(tmp, tmp, tmp[:,33], Int.(tmp[:,59]), true, arch)
end

@kernel function divide_half_kernel!(tmp, con)
    i = @index(Global, Linear)
    if con[i] == 1.0
        @inbounds tmp[i,4]  = tmp[i,5]  .* 0.45
        @inbounds tmp[i,5]  = tmp[i,5]  .* 0.45
        @inbounds tmp[i,6]  = tmp[i,6]  .* 0.45
        @inbounds tmp[i,7]  = tmp[i,7]  .* 0.5
        @inbounds tmp[i,8]  = tmp[i,8]  .* 0.5
        @inbounds tmp[i,9]  = tmp[i,9]  .* 0.5
        @inbounds tmp[i,10] = tmp[i,10] .* 0.5
        @inbounds tmp[i,11] = tmp[i,11] .+ 1.0
        @inbounds tmp[i,12] = 1.0
    end
    @inbounds tmp[i,33] = 0.0
end
function divide_half!(tmp, con, arch::Architecture)
    kernel! = divide_half_kernel!(device(arch), 256, (size(tmp,1),))
    event = kernel!(tmp, con)
    wait(device(arch), event)
    return nothing
end

