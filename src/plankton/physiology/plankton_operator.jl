##### copy active individuals to plank.dataₜ
@kernel function copyto_tmp_kernel!(plank, tmp, con, idx)
    i,j = @index(Global, NTuple)
    if con[i] == 1.0
        @inbounds tmp[idx[i],j] = copy(plank[i,j])
        @inbounds plank[i,j] = 0.0 # deactivate individual
    end
end
function copyto_tmp!(plank, tmp, con, idx::AbstractArray{Int64,1}, arch::Architecture)
    kernel! = copyto_tmp_kernel!(device(arch), (16,16), (size(plank,1),62))
    event = kernel!(plank, tmp, con, idx)
    wait(device(arch), event)
    return nothing
end

##### grazing and grazing loss
function grazing!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank[:,62] .= 0.0
    plank[:,62] .= cumsum(plank[:,31])
    ##### copy grazed individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank[:,31], Int.(plank[:,62]), arch)
    ##### calculate grazing loss
    calc_loss!(tmp, Int.(tmp[:,13:15]), arch, plk.DOC.data, plk.POC.data,
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.grazFracC, p.grazFracN, p.grazFracP, p.R_NC, p.R_PC)
end

##### mortality and mortality loss
function mortality!(plank, tmp, arch::Architecture, g::Grids, plk, p)
    ##### calculate index for timestepper.tmp
    plank[:,62] .= 0.0
    plank[:,62] .= cumsum(plank[:,32])
    ##### copy grazed individuals to timestepper.tmp
    copyto_tmp!(plank, tmp, plank[:,32], Int.(plank[:,62]), arch)
    ##### calculate grazing loss
    calc_loss!(tmp, Int.(tmp[:,13:15]), arch, plk.DOC.data, plk.POC.data,
               plk.DON.data, plk.PON.data, plk.DOP.data, plk.POP.data,
               g, p.mortFracC, p.mortFracN, p.mortFracP, p.R_NC, p.R_PC)
end

##### cell division
function divide_copy!(plank, tmp, arch::Architecture, plank_num::Int64)
    ##### calculate index for timestepper.tmp
    plank[:,62] .= 0.0
    plank[:,62] .= cumsum(plank[:,33])
    plank[:,62] .= plank[:,62] .+ plank_num
    copyto_tmp!(plank, tmp, plank[:,33], Int.(plank[:,62]), arch)
end

@kernel function divide_half_kernel!(tmp, con)
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
        @inbounds tmp[i,33] = 0.0
    end
end
function divide_half!(tmp, con, arch::Architecture)
    kernel! = divide_half_kernel!(device(arch), 256, (size(tmp,1),))
    event = kernel!(tmp, con)
    wait(device(arch), event)
    return nothing
end




# function divide!(plank, plank_num)
#     dvid_num = floor(Int64, sum(plank[:,33]))
#     if dvid_num > size(plank,1) - plank_num
#         throw(ArgumentError("INDIVIDUAL: individual number exceeded"))
#     end
#     if dvid_num > 0
#         plank[plank_num+1:plank_num+dvid_num, :] .= plank[Bool.(plank[:,33]), :]
#         plank[Bool.(plank[:,33]), 4]  .= plank[Bool.(plank[:,33]), 5]  .* 0.45
#         plank[Bool.(plank[:,33]), 5]  .= plank[Bool.(plank[:,33]), 5]  .* 0.45
#         plank[Bool.(plank[:,33]), 6]  .= plank[Bool.(plank[:,33]), 6]  .* 0.45
#         plank[Bool.(plank[:,33]), 7]  .= plank[Bool.(plank[:,33]), 7]  .* 0.5
#         plank[Bool.(plank[:,33]), 8]  .= plank[Bool.(plank[:,33]), 8]  .* 0.5
#         plank[Bool.(plank[:,33]), 9]  .= plank[Bool.(plank[:,33]), 9]  .* 0.5
#         plank[Bool.(plank[:,33]), 10] .= plank[Bool.(plank[:,33]), 10] .* 0.5
#         plank[Bool.(plank[:,33]), 11] .= plank[Bool.(plank[:,33]), 11] .+ 1.0
#         plank[Bool.(plank[:,33]), 12] .= 1.0
#     else
#         nothing
#     end
# end
