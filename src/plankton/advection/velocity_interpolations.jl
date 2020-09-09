##### find indices (halo points excluded)
@inline find_xF_ind(x, g::Grids) = x ≥ 0.0 ? x ÷ g.Δx + 1 : x ÷ g.Δx

@inline find_yF_ind(y, g::Grids) = y ≥ 0.0 ? y ÷ g.Δy + 1 : y ÷ g.Δy

@inline find_zF_ind(z, g::Grids) = (g.Nz * g.Δz + z) ≥ 0.0 ?
    (g.Nz * g.Δz + z) ÷ g.Δz + 1 : (g.Nz * g.Δz + z) ÷ g.Δz

@inline find_xC_ind(x, g::Grids) = (x - g.Δx * 0.5) ≥ 0.0 ?
    (x - g.Δx * 0.5) ÷ g.Δx + 1 : (x - g.Δx * 0.5) ÷ g.Δx

@inline find_yC_ind(y, g::Grids) = (y - g.Δy * 0.5) ≥ 0.0 ?
    (y - g.Δy * 0.5) ÷ g.Δy + 1 : (y - g.Δy * 0.5) ÷ g.Δy

@inline find_zC_ind(z, g::Grids) = (g.Nz * g.Δz + z - g.Δz * 0.5) ≥ 0.0 ?
    (g.Nz * g.Δz + z - g.Δz * 0.5 ÷ g.Δz + 1) : (g.Nz * g.Δz + z - g.Δz * 0.5) ÷ g.Δz

##### multi-linear interpolation
function trilinear_itpl(x, y, z, x₀, y₀, z₀, a, g::Grids)
    x₀ = Int(x₀) + g.Hx
    y₀ = Int(y₀) + g.Hy
    z₀ = Int(z₀) + g.Hz
    xᵈ = x - g.xF[x₀]
    yᵈ = y - g.yF[y₀]
    zᵈ = z - g.zF[z₀]
    vel_00 = a[x₀, y₀,     z₀    ] * (1 - xᵈ) + a[x₀ + 1, y₀,     z₀    ] * xᵈ
    vel_01 = a[x₀, y₀,     z₀ + 1] * (1 - xᵈ) + a[x₀ + 1, y₀,     z₀ + 1] * xᵈ
    vel_10 = a[x₀, y₀ + 1, z₀    ] * (1 - xᵈ) + a[x₀ + 1, y₀ + 1, z₀    ] * xᵈ
    vel_11 = a[x₀, y₀ + 1, z₀ + 1] * (1 - xᵈ) + a[x₀ + 1, y₀ + 1, z₀ + 1] * xᵈ
    vel_0 = vel_00 * (1 - yᵈ) + vel_10 * yᵈ
    vel_1 = vel_01 * (1 - yᵈ) + vel_11 * yᵈ
    vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ
    return vel
end
