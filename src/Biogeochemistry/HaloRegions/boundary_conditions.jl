#####
##### generate boundary conditions for each nutrient field
#####

# generate the default boundary conditions for one field
function default_bcs()
    bcs = BoundaryConditions(nothing, nothing, nothing, nothing, nothing, nothing)
    return bcs
end

"""
    set_bc!(model; tracer::Symbol, pos::Symbol, bc_value::Union{Number, AbstractArray})
Set the boundary condition of `tracer` on `pos` with `bc_value` of DataType `FT`.

Keyword Arguments
=================
- `tracer`: the tracer of which the boundary condition will be set.
- `pos`: the position of the bounday condition to be set, e.g., `:east`, `:top` etc.
- `bc_value`: the value that will be used to set the boundary condition.
"""
function set_bc!(model; tracer::Symbol, pos::Symbol, bc_value::Union{Number, AbstractArray})
    @assert tracer in tracer_names

    FT = model.FT
    bc_value_d = FT.(bc_value)
    if isa(bc_value_d, AbstractArray)
        bc_value_d = bc_value_d |> array_type(model.arch)
    end
    setproperty!(model.tracers[tracer].bc, pos, bc_value_d)
    return nothing
end

# get boundary condition at each grid point
@inline function getbc(bc::Union{Number, AbstractArray}, i, j, t)
    if typeof(bc) <: Number
        return bc
    elseif typeof(bc) <: AbstractArray{eltype(bc),2}
        return bc[i,j]
    elseif typeof(bc) <: AbstractArray{eltype(bc),3}
        return bc[i,j,t]
    end
end

# validate boundary conditions, check if the grid information is compatible with nutrient field
function validate_bc(bc::Union{Nothing, Number, AbstractArray}, bc_size, nΔT)
    if typeof(bc) <: AbstractArray{eltype(bc),2} 
        if size(bc) == bc_size
            return nothing
        else
            throw(ArgumentError("BC: grid mismatch, size(bc) must equal to $(bc_size) for a constant flux boundary condition."))
        end
    end
    if typeof(bc) <: AbstractArray{eltype(bc),3}
        if size(bc) == (bc_size..., nΔT)
            return nothing
        else
            throw(ArgumentError("BC: grid mismatch, size(bc) must equal to $((bc_size..., nΔT)) for a time-dependent flux boundary condition."))
        end
    end
end

function validate_bcs(nut, grid, nΔT)
    for name in keys(nut)
        bc_west    = nut[name].bc.west
        bc_east    = nut[name].bc.east
        bc_south   = nut[name].bc.south
        bc_north   = nut[name].bc.north
        bc_bottom  = nut[name].bc.bottom
        bc_top     = nut[name].bc.top

        bc_size_x = (grid.Ny, grid.Nz)
        bc_size_y = (grid.Nx, grid.Nz)
        bc_size_z = (grid.Nx, grid.Ny)
        
        validate_bc(bc_west, bc_size_x, nΔT)
        validate_bc(bc_east, bc_size_x, nΔT)
        validate_bc(bc_south, bc_size_y, nΔT)
        validate_bc(bc_north, bc_size_y, nΔT)
        validate_bc(bc_bottom, bc_size_z, nΔT)
        validate_bc(bc_top, bc_size_z, nΔT)
    end
end





