##### calculate diffusivities of each individual
function plankton_diffusion!(plank, rnd, κx, κy, κz, ΔT, num::Int64)
    plank[1:num,1] .= plank[1:num,1] .+ rnd[1:num,1] .* √(3*ΔT) .* (2*√(2*κx))
    plank[1:num,2] .= plank[1:num,2] .+ rnd[1:num,2] .* √(3*ΔT) .* (2*√(2*κy))
    plank[1:num,3] .= plank[1:num,3] .+ rnd[1:num,3] .* √(3*ΔT) .* (2*√(2*κz))

    return nothing
end

function gen_rand_adv!(rnd, arch)
    rand!(rng_type(arch), rnd)
    rnd .= rnd .* 2.0 .- 1.0
end

plankton_diffusion!(plank, rnd, κ, ΔT, num) = plankton_diffusion!(plank, rnd, κ, κ, κ, ΔT, num)
