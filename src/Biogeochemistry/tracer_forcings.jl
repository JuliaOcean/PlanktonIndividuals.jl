##### calculate simple remineralization of DOM and POM as well as a simple nitrification
function tracer_forcing!(F, tracer_temp, tracer, params, ΔT)
    for name in tracer_names
        @inbounds tracer_temp[name].data .= max.(0.0f0, tracer[name].data)
    end

    @inbounds F.DIC.data .= F.DIC.data .+ tracer_temp.DOC.data .* params["kDOC"] .* ΔT

    @inbounds F.DOC.data .= F.DOC.data .- tracer_temp.DOC.data .* params["kDOC"] .* ΔT .+
                                          tracer_temp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.POC.data .= F.POC.data .- tracer_temp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.NH4.data .= F.NH4.data .+ tracer_temp.DON.data .* params["kDON"] .* ΔT .-
                                          tracer_temp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.NO3.data .= F.NO3.data .+ tracer_temp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.DON.data .= F.DON.data .- tracer_temp.DON.data .* params["kDON"] .* ΔT .+
                                          tracer_temp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PON.data .= F.PON.data .- tracer_temp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PO4.data .= F.PO4.data .+ tracer_temp.DOP.data .* params["kDOP"] .* ΔT

    @inbounds F.DOP.data .= F.DOP.data .- tracer_temp.DOP.data .* params["kDOP"] .* ΔT .+
                                          tracer_temp.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.POP.data .= F.POP.data .- tracer_temp.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.FeT.data .= F.FeT.data .+ tracer_temp.DOFe.data .* params["kDOFe"] .* ΔT

    @inbounds F.DOFe.data .= F.DOFe.data .- tracer_temp.DOFe.data .* params["kDOFe"] .* ΔT .+
                                          tracer_temp.POFe.data .* params["kPOFe"] .* ΔT

    @inbounds F.POFe.data .= F.POFe.data .- tracer_temp.POFe.data .* params["kPOFe"] .* ΔT

    @inbounds F.CHO.data .= F.CHO.data .- tracer_temp.CHO.data .* params["kCHO"] .* ΔT

    return nothing
end

