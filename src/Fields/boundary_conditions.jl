#####
##### generate boundary conditions for each nutrient field
#####

# struct to combine the left and right boundary conditions in one direction
mutable struct BoundaryConditions
    west::Union{Nothing, Number, AbstractArray}
    east::Union{Nothing, Number, AbstractArray}
    north::Union{Nothing, Number, AbstractArray}
    south::Union{Nothing, Number, AbstractArray}
    top::Union{Nothing, Number, AbstractArray}
    bottom::Union{Nothing, Number, AbstractArray}
end

# generate the default boundary conditions for one field
function default_bcs()
    bcs = BoundaryConditions(nothing, nothing, nothing, nothing, nothing, nothing)
    return bcs
end

"""
    set_bc!(model, tracer::Symbol, pos::Symbol, bc_value::Union{Number, AbstractArray})
Set the boundary condition of `tracer` on `pos` with `bc_value`.
"""
function set_bc!(model, tracer::Symbol, pos::Symbol, bc_value::Union{Number, AbstractArray})
    @assert tracer in nut_names

    bc_value_d = bc_value
    if isa(bc_value, AbstractArray)
        bc_value_d = vc_value |> array_type(model.arch)
    end
    setproperty!(model.nutrients[tracer].bc, pos, bc_value_d)
end

# get boundary condition at each grid point
@inline getbc(bc::Number, i, j, t) = bc
@inline getbc(bc::AbstractArray{Float64,2}, i, j, t) = bc[i,j]
@inline getbc(bc::AbstractArray{Float64,3}, i, j, t) = bc[i,j,t]

# validate boundary conditions, check if the grid information is compatible with nutrient field
function validate_bc(bc, bc_size, nΔT)
    if typeof(bc) <: AbstractArray{Float64,2} 
        if size(bc) == bc_size
            return nothing
        else
            throw(ArgumentError("BC west: grid mismatch, size(bc) must equal to $(bc_size) for a constant flux boundary condition."))
        end
    end
    if typeof(bc) <: AbstractArray{Float64,3}
        if size(bc) == (bc_size..., nΔT)
            return nothing
        else
            throw(ArgumentError("BC west: grid mismatch, size(bc) must equal to $((bc_size..., nΔT)) for a time-dependent flux boundary condition."))
        end
    end
end

function validate_bcs(nut, grid, nΔT)
    for name in nut_names
        bc_west    = nut[name].bc.west
        bc_east    = nut[name].bc.east
        bc_south   = nut[name].bc.south
        bc_north   = nut[name].bc.north
        bc_bottom  = nut[name].bc.bottom
        bc_top     = nut[name].bc.top

        bc_size_x = (grid.Ny+2*grid.Hy, grid.Nz+2*grid.Hz)
        bc_size_y = (grid.Nx+2*grid.Hx, grid.Nz+2*grid.Hz)
        bc_size_z = (grid.Nx+2*grid.Hx, grid.Ny+2*grid.Hy)
        
        validate_bc(bc_west, bc_size_x, nΔT)
        validate_bc(bc_east, bc_size_x, nΔT)
        validate_bc(bc_south, bc_size_y, nΔT)
        validate_bc(bc_north, bc_size_y, nΔT)
        validate_bc(bc_bottom, bc_size_z, nΔT)
        validate_bc(bc_top, bc_size_z, nΔT)
    end
end





