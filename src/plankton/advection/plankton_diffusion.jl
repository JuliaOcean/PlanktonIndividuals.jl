##### calculate diffusivities of each individual
function plankton_diffusion!(plank, rnd, κx, κy, κz, ΔT)
    plank.x .= plank.x .+ rnd.x .* √(3*ΔT) .* (2*√(2*κx)) .* plank.ac
    plank.y .= plank.y .+ rnd.y .* √(3*ΔT) .* (2*√(2*κy)) .* plank.ac
    plank.z .= plank.z .+ rnd.z .* √(3*ΔT) .* (2*√(2*κz)) .* plank.ac

    return nothing
end

plankton_diffusion!(plank, rnd, κ, ΔT) = plankton_diffusion!(plank, rnd, κ, κ, κ, ΔT)

function gen_rand_adv!(rnd, arch)
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)
    rnd.x .= rnd.x .* 2.0 .- 1.0
    rnd.y .= rnd.y .* 2.0 .- 1.0
    rnd.z .= rnd.z .* 2.0 .- 1.0
end
