module Fields

export Field
export interior, zero_fields!
export fill_halo_nut!, fill_halo_Gcs!, fill_halo_vel!
export fill_halo_u!, fill_halo_v!, fill_halo_w!
export apply_bcs!
export default_bcs, getbc, validate_bcs, LFBoundaryConditions
export nut_names

using KernelAbstractions
using CUDA

using PlanktonIndividuals.Grids
using PlanktonIndividuals.Architectures: device, Architecture, array_type


struct Field
    data::AbstractArray{Float64,3}
    bc::NamedTuple
end

function Field(arch::Architecture, g::AbstractGrid; bcs = default_bcs())
    total_size = (g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2)
    data = zeros(total_size) |> array_type(arch)
    return Field(data,bcs)
end

@inline interior(c, g) = c[g.Hx+1:g.Hx+g.Nx, g.Hy+1:g.Hy+g.Ny, g.Hz+1:g.Hz+g.Nz]

function zero_fields!(a)
    for i in 1:length(a)
        @inbounds a[i].data .= 0.0
    end
end

const nut_names=(:DIC,:NH4,:NO3,:PO4,:DOC,:DON,:DOP,:POC,:PON,:POP)

include("halo_regions.jl")
include("boundary_conditions.jl")
include("apply_bcs.jl")

end