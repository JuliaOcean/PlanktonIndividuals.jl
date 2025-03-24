##### calculate simple remineralization of DOM and POM as well as a simple nitrification
function nut_forcing!(F, nut_temp, nut, params, ΔT)
    for name in nut_names
        @inbounds nut_temp[name].data .= max.(0.0f0, nut[name].data)
    end

    @inbounds F.DIC.data .= F.DIC.data .+ nut_temp.DOC.data .* params["kDOC"] .* ΔT

    @inbounds F.DOC.data .= F.DOC.data .- nut_temp.DOC.data .* params["kDOC"] .* ΔT .+
                                          nut_temp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.POC.data .= F.POC.data .- nut_temp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.NH4.data .= F.NH4.data .+ nut_temp.DON.data .* params["kDON"] .* ΔT .-
                                          nut_temp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.NO3.data .= F.NO3.data .+ nut_temp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.DON.data .= F.DON.data .- nut_temp.DON.data .* params["kDON"] .* ΔT .+
                                          nut_temp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PON.data .= F.PON.data .- nut_temp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PO4.data .= F.PO4.data .+ nut_temp.DOP.data .* params["kDOP"] .* ΔT

    @inbounds F.DOP.data .= F.DOP.data .- nut_temp.DOP.data .* params["kDOP"] .* ΔT .+
                                          nut_temp.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.POP.data .= F.POP.data .- nut_temp.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.FeT.data .= F.FeT.data .+ nut_temp.DOFe.data .* params["kDOFe"] .* ΔT

    @inbounds F.DOFe.data .= F.DOFe.data .- nut_temp.DOFe.data .* params["kDOFe"] .* ΔT .+
                                          nut_temp.POFe.data .* params["kPOFe"] .* ΔT

    @inbounds F.POFe.data .= F.POFe.data .- nut_temp.POFe.data .* params["kPOFe"] .* ΔT

    @inbounds F.CHO.data .= F.CHO.data .- nut_temp.CHO.data .* params["kCHO"] .* ΔT

    return nothing
end

