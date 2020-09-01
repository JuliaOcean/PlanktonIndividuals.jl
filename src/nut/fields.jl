struct Field
    data::AbstractArray{Float64,3}
    grid::grids
end

function Field(::CPUs,grid)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    data = zeros(total_size)
    return Field(data,grid)
end

function Field(::GPUs,grid)
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2+1)
    data = zeros(total_size) |> CuArray
    return Field(data,grid)
end

const nut_names=(:DIC,:NH4,:NO3,:PO4,:DOC,:DON,:DOP,:POC,:PON,:POP)
