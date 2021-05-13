struct Field
    data::AbstractArray{Float64,3}
    bc::NamedTuple
end

function Field(arch::Architecture, g::AbstractGrid)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    data = zeros(total_size) |> array_type(arch)
    bcs = default_bcs()
    return Field(data,bcs)
end

const nut_names=(:DIC,:NH4,:NO3,:PO4,:DOC,:DON,:DOP,:POC,:PON,:POP)

@inline interior(c, g) = c[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz]
