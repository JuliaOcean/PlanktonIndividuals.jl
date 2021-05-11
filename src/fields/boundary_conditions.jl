#####
##### generate boundary conditions for each nutrient field
#####

# struct to combine the left and right boundary conditions in one direction
mutable struct LFBoundaryConditions
    left::Union{Nothing, Number, AbstractArray}
    right::Union{Nothing, Number, AbstractArray}
end

# generate the default boundary conditions for one field
function default_bcs()
    x = LFBoundaryConditions(nothing, nothing)
    y = LFBoundaryConditions(nothing, nothing)
    z = LFBoundaryConditions(nothing, nothing)

    bcs = (x = x, y = y, z = z)

    return bcs
end

# get boundary condition at each grid point
@inline getbc(bc::Number, i, j, t) = bc
@inline getbc(bc::AbstractArray{Float64,2}, i, j, t) = bc[i,j]
@inline getbc(bc::AbstractArray{Float64,3}, i, j, t) = bc[i,j,t]