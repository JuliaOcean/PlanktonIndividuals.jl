##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ 0.0 ? Int(x ÷ g.Δx + 1) : Int(x ÷ g.Δx)

@inline find_yF_ind(y, g::Grids) = y ≥ 0.0 ? Int(y ÷ g.Δy + 1) : Int(y ÷ g.Δy)

@inline find_zF_ind(z, g::Grids) = (g.Nz * g.Δz + z) ≥ 0.0 ?
    Int((g.Nz * g.Δz + z) ÷ g.Δz + 1) : Int((g.Nz * g.Δz + z) ÷ g.Δz)

@inline find_xC_ind(x, g::Grids) = (x - g.Δx * 0.5) ≥ 0.0 ?
    Int((x - g.Δx * 0.5) ÷ g.Δx + 1) : Int((x - g.Δx * 0.5) ÷ g.Δx)

@inline find_yC_ind(y, g::Grids) = (y - g.Δy * 0.5) ≥ 0.0 ?
    Int((y - g.Δy * 0.5) ÷ g.Δy + 1) : Int((y - g.Δy * 0.5) ÷ g.Δy)

@inline find_zC_ind(z, g::Grids) = (g.Nz * g.Δz + z - g.Δz * 0.5) ≥ 0.0 ?
    Int((g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz + 1) : Int((g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz)

##### multi-linear interpolation
function trilinear_itpl(coords, x₀, y₀, z₀, a)
    xᵈ = coords[1] - x₀
    yᵈ = coords[2] - y₀
    zᵈ = coords[3] - z₀
    vel_00 = a[x₀, y₀,     z₀    ] * (1 - yᵈ) + a[x₀ + 1, y₀,     z₀    ] * yᵈ
    vel_01 = a[x₀, y₀,     z₀ + 1] * (1 - yᵈ) + a[x₀ + 1, y₀,     z₀ + 1] * yᵈ
    vel_10 = a[x₀, y₀ + 1, z₀    ] * (1 - yᵈ) + a[x₀ + 1, y₀ + 1, z₀    ] * yᵈ
    vel_11 = a[x₀, y₀ + 1, z₀ + 1] * (1 - yᵈ) + a[x₀ + 1, y₀ + 1, z₀ + 1] * yᵈ
    vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
    vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
    vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ
    return vel
end
