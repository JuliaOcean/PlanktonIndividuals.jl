##### calculate simple remineralization of DOM and POM as well as a simple nitrification
function nut_forcing!(F, MD1, nut, params, ΔT)
    for name in nut_names
        @inbounds MD1[name].data .= max.(0.0, nut[name].data)
    end

    @inbounds F.DIC.data .= F.DIC.data .+ MD1.DOC.data .* params["kDOC"] .* ΔT

    @inbounds F.DOC.data .= F.DOC.data .- MD1.DOC.data .* params["kDOC"] .* ΔT .+
                                          MD1.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.POC.data .= F.POC.data .- MD1.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.NH4.data .= F.NH4.data .+ MD1.DON.data .* params["kDON"] .* ΔT .-
                                          MD1.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.NO3.data .= F.NO3.data .+ MD1.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.DON.data .= F.DON.data .- MD1.DON.data .* params["kDON"] .* ΔT .+
                                          MD1.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PON.data .= F.PON.data .- MD1.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PO4.data .= F.PO4.data .+ MD1.DOP.data .* params["kDOP"] .* ΔT

    @inbounds F.DOP.data .= F.DOP.data .- MD1.DOP.data .* params["kDOP"] .* ΔT .+
                                          MD1.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.POP.data .= F.POP.data .- MD1.POP.data .* params["kPOP"] .* ΔT

    return nothing
end

