# for offline only, include time series
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

struct grid
    xC::Array
    yC::Array
    zC::Array
    xF::Array
    yF::Array
    zF::Array
    Δx::Array # converted to m
    Δy::Array # converted to m
    Δz::Array # converted to m
    Nx::Int
    Ny::Int
    Nz::Int
end

mutable struct nutreint_fields
    DIC::Array{Float64,3}
    DIN::Array{Float64,3}
    DOC::Array{Float64,3}
    DON::Array{Float64,3}
    POC::Array{Float64,3}
    PON::Array{Float64,3}
end
