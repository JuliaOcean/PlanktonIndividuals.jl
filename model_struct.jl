#for offline only, include time series
struct velocity_fields 
    u::Array{Float64,4}
    v::Array{Float64,4}
    w::Array{Float64,4}
end

mutable struct velocity
    u::Array{Float64,3}
    v::Array{Float64,3}
    w::Array{Float64,3}
end

struct grids
    xC::Array
    yC::Array
    zC::Array
    xF::Array
    yF::Array
    zF::Array
    Δx::Array # unit: degree
    Δy::Array # unit: degree
    Δz::Array # unit: degree
    Nx::Int
    Ny::Int
    Nz::Int
end

mutable struct nutrient_fields # for tendencies, forcings and consumptions
    DIC::Array{Float64,3}
    DIN::Array{Float64,3}
    DOC::Array{Float64,3}
    DON::Array{Float64,3}
    POC::Array{Float64,3}
    PON::Array{Float64,3}
end

struct rem # parameters for nutrient remineralization
    DOC::Float64
    DON::Float64
    POC::Float64
    PON::Float64
end
