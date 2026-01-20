abstract type AbstractTankController end

# Mutable struct to hold the cache (because optimizer_cache changes)
mutable struct MPCController <: AbstractTankController
    params::TankParameters
    horizon::Int
    references::Vector{Float64}
    optimizer_cache::Union{MPCModelRefs,Nothing}
    last_u::Dict{Symbol,Any}
    lambda_switch::Float64
    settings::InputSettings
end

function MPCController(;
    params=TankParameters(),
    horizon=10,
    references=[0.5, 0.5, 0.5],
    lambda_switch=0.1,
    settings=InputSettings()
)
    return MPCController(params, horizon, references, nothing, Dict{Symbol,Any}(), lambda_switch, settings)
end


function compute_control(c::MPCController, current_h::Vector, t::Float64)
    # Rebuild model every step to avoid bridge issues
    c.optimizer_cache = build_mpc_model(c.params, c.horizon, c.lambda_switch)

    # Determine fixed inputs from settings
    # All inputs - enabled_inputs = fixed inputs
    # fixed value comes from nominal_values
    all_inputs = [:Q1, :Q2, :V1, :V2, :V13, :V23, :VL1, :VL2, :VL3]
    fixed_u_map = Dict{Symbol,Any}()

    for inp in all_inputs
        if !(inp in c.settings.enabled_inputs)
            # Input is NOT enabled -> FIX IT
            val = get(c.settings.nominal_values, inp, 0.0)
            fixed_u_map[inp] = val
        end
    end

    update_mpc_model!(c.optimizer_cache, current_h, c.references, c.last_u, fixed_u_map)

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
