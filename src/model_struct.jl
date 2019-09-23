#for offline only, include time series
mutable struct velocity
    u::Array{Float64,3}
    v::Array{Float64,3}
    w::Array{Float64,3}
end

struct grids
    xC::Array{Float32,2}
    yC::Array{Float32,2}
    zC::Array{Float32,1}
    xF::Array{Float32,2}
    yF::Array{Float32,2}
    zF::Array{Float32,1}
    Δx::Array{Float32,2} # unit: degree
    Δy::Array{Float32,2} # unit: degree
    Lx::Array{Float32,2} # unit: meter
    Ly::Array{Float32,2} # unit: meter
    Lz::Array{Float32,1} # unit: meter
    dxC::Array{Float32,2} # unit: meter, distance from center to center
    dyC::Array{Float32,2} # unit: meter, distance from center to center
    dzC::Array{Float32,1} # unit: meter, distance from center to center
    Ax::Array{Float32,3} # unit: m²
    Ay::Array{Float32,3} # unit: m²
    Az::Array{Float32,2} # unit: m²
    V ::Array{Float32,3} # unit: m³
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

struct remineralization # parameters for nutrient remineralization
    DOC::Float64
    DON::Float64
    POC::Float64
    PON::Float64
end
