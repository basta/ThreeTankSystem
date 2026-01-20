abstract type AbstractTankController end

# Mutable struct to hold the cache (because optimizer_cache changes)
mutable struct MPCController <: AbstractTankController
    params::TankParameters
    horizon::Int
    references::Vector{Float64}
    optimizer_cache::Union{MPCModelRefs,Nothing}
    last_u::Dict{Symbol,Any}
end

function MPCController(;
    params=TankParameters(),
    horizon=10,
    references=[0.5, 0.5, 0.5]
)
    return MPCController(params, horizon, references, nothing, Dict{Symbol,Any}())
end


function compute_control(c::MPCController, current_h::Vector, t::Float64; fixed_u=Dict())
    if isnothing(c.optimizer_cache)
        c.optimizer_cache = build_mpc_model(c.params, c.horizon)
    end

    update_mpc_model!(c.optimizer_cache, current_h, c.references, c.last_u, fixed_u)

    optimize!(c.optimizer_cache.model)
    u_opt = extract_solution(c.optimizer_cache)

    c.last_u = Dict(
        :Q1 => u_opt.Q1,
        :Q2 => u_opt.Q2,
        :V1 => u_opt.V1,
        :V2 => u_opt.V2,
        :V13 => u_opt.V13,
        :V23 => u_opt.V23,
        :VL1 => u_opt.VL1,
        :VL2 => u_opt.VL2,
        :VL3 => u_opt.VL3
    )

    return u_opt
end
