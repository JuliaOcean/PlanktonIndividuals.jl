module Fields

export Field
export interior, zero_fields!
export fill_halo_nut!, fill_halo_Gcs!, fill_halo_vel!
export fill_halo_u!, fill_halo_v!, fill_halo_w!
export apply_bcs!, getbc
export default_bcs
export set_bc!, validate_bcs
export nut_names

using KernelAbstractions
using CUDA

using PlanktonIndividuals.Grids
using PlanktonIndividuals.Architectures: device, Architecture, array_type

include("halo_regions.jl")
include("boundary_conditions.jl")
include("apply_bcs.jl")

struct Field{FT}
    data::AbstractArray{FT,3}
    bc::BoundaryConditions
end

"""
    Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
Construct a `Field` on `grid` with data and boundary conditions on architecture `arch`
with DataType `FT`.
"""
function Field(arch::Architecture, grid::AbstractGrid, FT::DataType; bcs = default_bcs())
    total_size = (grid.Nx+grid.Hx*2, grid.Ny+grid.Hy*2, grid.Nz+grid.Hz*2)
    data = zeros(FT, total_size) |> array_type(arch)
    return Field{FT}(data,bcs)
end

@inline interior(c, grid) = c[grid.Hx+1:grid.Hx+grid.Nx, grid.Hy+1:grid.Hy+grid.Ny, grid.Hz+1:grid.Hz+grid.Nz]

function zero_fields!(a)
    for i in 1:length(a)
        @inbounds a[i].data .= 0.0f0
    end
end

const nut_names=(:DIC,:NH4,:NO3,:PO4,:DOC,:DON,:DOP,:POC,:PON,:POP)

end