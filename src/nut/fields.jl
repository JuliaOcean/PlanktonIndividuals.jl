struct Field
    data::AbstractArray{Float64,3}
    grid::Grids
end

function Field(::CPUs, g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    data = zeros(total_size)
    return Field(data,g)
end

function Field(::GPUs, g)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    data = zeros(total_size) |> CuArray
    return Field(data,g)
end

const nut_names=(:DIC,:NH4,:NO3,:PO4,:DOC,:DON,:DOP,:POC,:PON,:POP)

@inline interior(c, g) = c[g.Hx+1:g.Hx+g.Nz, g.Hy+1:g.Hy+g.Hy, g.Hz+1:g.Hz+g.Nz]
