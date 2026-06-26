##### materials exchange between cells with a colony
##### O₂, qP, qFe, ATP, NADPH, CH, NH₃
@inline function calc_qP_exchange(qP1, qP2, sp1_p)
    ΔqP = qP1 - qP2
    TqP = -ΔqP * sp1_p.k_TqP * (1 - sp1_p.ϵ_TqP) * sp1_p.Nsuper
    return TqP
end
@inline function calc_CH_exchange(CH1, CH2, sp1_p)
    ΔCH = CH1 - CH2
    TCH = -ΔCH * sp1_p.k_TCH * (1 - sp1_p.ϵ_TCH) * sp1_p.Nsuper
    return TCH
end
@inline function calc_qO2_exchange(qO21, qO22, sp1_p)
    ΔqO2 = qO21 - qO22
    TqO2 = -ΔqO2 * sp1_p.k_TqO2 * (1 - sp1_p.ϵ_TqO2) * sp1_p.Nsuper
    return TqO2
end
@inline function calc_qFe_exchange(qFe1, qFe2, sp1_p)
    ΔqFe = qFe1 - qFe2
    TqFe = -ΔqFe * sp1_p.k_TqFe * (1 - sp1_p.ϵ_TqFe) * sp1_p.Nsuper
    return TqFe
end
@inline function calc_NH4_exchange(qNH41, qNH42, sp1_p)
    ΔNH3 = qNH41 - qNH42
    TNH3 = -ΔNH3 * sp1_p.k_TqNH4 * (1 - sp1_p.ϵ_TqNH4) * sp1_p.Nsuper
    return TNH3
end


@kernel function calc_material_exchange_kernel!(sp1, sp2, sp1_p)
    i = @index(Global)
    @inbounds sp1.TqP[i]   = calc_qP_exchange(sp1.qP[i], sp2.qP[i], sp1_p) * sp1.ac[i]
    @inbounds sp1.TCH[i]   = calc_CH_exchange(sp1.CH[i], sp2.CH[i], sp1_p) * sp1.ac[i]
    @inbounds sp1.TqO2[i]  = calc_qO2_exchange(sp1.qO2[i], sp2.qO2[i], sp1_p) * sp1.ac[i]
    @inbounds sp1.TqFe[i]  = calc_qFe_exchange(sp1.qFe[i], sp2.qFe[i], sp1_p) * sp1.ac[i]
    @inbounds sp1.TqNH4[i] = calc_NH4_exchange(sp1.qNH4[i], sp2.qNH4[i], sp1_p) * sp1.ac[i]
    @inbounds sp2.TqP[i]   = -sp1.TqP[i]
    @inbounds sp2.TCH[i]   = -sp1.TCH[i]
    @inbounds sp2.TqO2[i]  = -sp1.TqO2[i]
    @inbounds sp2.TqFe[i]  = -sp1.TqFe[i]
    @inbounds sp2.TqNH4[i] = -sp1.TqNH4[i]
end
function calc_material_exchange!(sp1, sp2, sp1_p, arch::Architecture)
    kernel! = calc_material_exchange_kernel!(device(arch), 256, (size(sp1.ac,1)))
    kernel!(sp1, sp2, sp1_p)
    return nothing
end

@kernel function update_material_exchange_kernel!(sp1, sp2, ΔT)
    i = @index(Global)
    @inbounds sp1.qP[i]    += ΔT * sp1.TqP[i]
    @inbounds sp2.qP[i]    += ΔT * sp2.TqP[i]
    @inbounds sp1.CH[i]    += ΔT * sp1.TCH[i]
    @inbounds sp2.CH[i]    += ΔT * sp2.TCH[i]
    @inbounds sp1.qO2[i]   += ΔT * sp1.TqO2[i]
    @inbounds sp2.qO2[i]   += ΔT * sp2.TqO2[i]
    @inbounds sp1.qNH4[i]  += ΔT * sp1.TqNH4[i]
    @inbounds sp2.qNH4[i]  += ΔT * sp2.TqNH4[i]
    @inbounds sp1.qFe[i]   += ΔT * sp1.TqFe[i]
    @inbounds sp2.qFe[i]   += ΔT * sp2.TqFe[i]
end
function update_material_exchange!(sp1, sp2, ΔT, arch::Architecture)
    kernel! = update_material_exchange_kernel!(device(arch), 256, (size(sp1.ac,1)))
    kernel!(sp1, sp2, ΔT)
    return nothing
end