"""
    function nut_update(model, velᵇ, consume, ΔT, Dim)
Update nutrient fields to next time step with source term and consumption by phytoplankton ('consume')
"""
function nut_update(model, velᵇ, consume, ΔT, Dim=3)
    nutrients = model.nutrients
    g = model.grid
    κh = model.params["κh"]
    κv = model.params["κv"]
    # compute Forcings (remineralization)
    DIC = copy(nutrients.DIC);DIN = copy(nutrients.DIN);
    DOC = copy(nutrients.DOC);DON = copy(nutrients.DON);
    POC = copy(nutrients.POC);PON = copy(nutrients.PON);
    DIC[DIC .< 0.0] .= 0.0;DIN[DIN .< 0.0] .= 0.0;
    DOC[DOC .< 0.0] .= 0.0;DON[DON .< 0.0] .= 0.0;
    POC[POC .< 0.0] .= 0.0;PON[PON .< 0.0] .= 0.0;
    F = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                        zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                        zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    # compute remineralization of organic nutrients
    F.DIC .= F.DIC .+ DOC .* model.params["kDOC"]
    F.DIN .= F.DIN .+ DON .* model.params["kDON"]
    F.DOC .= F.DOC .- DOC .* model.params["kDOC"] .+ POC .* model.params["kPOC"]
    F.DON .= F.DON .- DON .* model.params["kDON"] .+ PON .* model.params["kPON"]
    F.POC .= F.POC .- POC .* model.params["kPOC"]
    F.PON .= F.PON .- PON .* model.params["kPON"]

    # Compute nutrient advection using DST3FL scheme 'MultiDim_adv(g, DIC, vel, ΔT)
    # Compute nutrient diffusion 'κ∇²(g, DIC, κh, κv, i, j, k)
    gtr = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                          zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    if Dim == 0 
        nothing
    else
        adv_DIC = MultiDim_adv(g, DIC, velᵇ, ΔT);
        adv_DIN = MultiDim_adv(g, DIN, velᵇ, ΔT);
        adv_DOC = MultiDim_adv(g, DOC, velᵇ, ΔT);
        adv_DON = MultiDim_adv(g, DON, velᵇ, ΔT);
        adv_POC = MultiDim_adv(g, POC, velᵇ, ΔT);
        adv_PON = MultiDim_adv(g, PON, velᵇ, ΔT);
        for k in 1:g.Nz
            for j in 1:g.Ny
                for i in 1:g.Nx
                    gtr.DIC[i, j, k] = adv_DIC[i, j, k] + κ∇²(g, DIC, model.params["κh"], model.params["κh"], i, j, k) + F.DIC[i, j, k]
                    gtr.DIN[i, j, k] = adv_DIN[i, j, k] + κ∇²(g, DIN, model.params["κh"], model.params["κh"], i, j, k) + F.DIN[i, j, k]
                    gtr.DOC[i, j, k] = adv_DOC[i, j, k] + κ∇²(g, DOC, model.params["κh"], model.params["κh"], i, j, k) + F.DOC[i, j, k]
                    gtr.DON[i, j, k] = adv_DON[i, j, k] + κ∇²(g, DON, model.params["κh"], model.params["κh"], i, j, k) + F.DON[i, j, k]
                    gtr.POC[i, j, k] = adv_POC[i, j, k] + κ∇²(g, POC, model.params["κh"], model.params["κh"], i, j, k) + F.POC[i, j, k]
                    gtr.PON[i, j, k] = adv_PON[i, j, k] + κ∇²(g, PON, model.params["κh"], model.params["κh"], i, j, k) + F.PON[i, j, k]
                end
            end
        end
    end
    # store nutrients of the former time step
    nutₜ = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                           zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz),
                           zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
        nutₜ.DIC .= nutrients.DIC .+ gtr.DIC .* ΔT .+ consume.DIC ./ g.V
        nutₜ.DIN .= nutrients.DIN .+ gtr.DIC .* ΔT .+ consume.DIN ./ g.V
        nutₜ.DOC .= nutrients.DOC .+ gtr.DOC .* ΔT .+ consume.DOC ./ g.V
        nutₜ.DON .= nutrients.DON .+ gtr.DON .* ΔT .+ consume.DON ./ g.V
        nutₜ.POC .= nutrients.POC .+ gtr.POC .* ΔT .+ consume.POC ./ g.V
        nutₜ.PON .= nutrients.PON .+ gtr.PON .* ΔT .+ consume.PON ./ g.V
    return nutₜ, gtr
end

