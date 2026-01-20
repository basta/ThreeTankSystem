using OrdinaryDiffEq

mutable struct Simulation
    plant::ODESystem
    controller::AbstractTankController
    current_state::Vector{Float64}
    history::Dict{Symbol,Any}
    p_static::Dict{Any,Any}
    p_nominal::Dict{Any,Any}
end

function Simulation(; plant, controller, h0)
    p = controller.params

    # Helper to robustly set parameters if they exist in the system
    function try_set!(d, sys, name, val)
        try
            # getproperty throws if name is not a variable/parameter in sys
            sym = getproperty(sys, name)
            d[sym] = val
        catch
        end
    end

    p_static = Dict{Any,Any}()

    # Map common parameters
    try_set!(p_static, plant, :A, p.A)
    try_set!(p_static, plant, :az, p.az)
    try_set!(p_static, plant, :g, p.g)
    try_set!(p_static, plant, :hv, p.hv)
    try_set!(p_static, plant, :hmax, p.hmax)
    try_set!(p_static, plant, :Qimax, p.Qimax)

    try_set!(p_static, plant, :SL1, p.SL1)
    try_set!(p_static, plant, :SL2, p.SL2)
    try_set!(p_static, plant, :SL3, p.SL3)
    try_set!(p_static, plant, :SL13, p.SL13)
    try_set!(p_static, plant, :SL23, p.SL23)
    try_set!(p_static, plant, :S1, p.S1)
    try_set!(p_static, plant, :S2, p.S2)

    # Nominal inputs (default closed)
    p_nominal = Dict{Any,Any}()
    try_set!(p_nominal, plant, :Q1, 0.0)
    try_set!(p_nominal, plant, :Q2, 0.0)
    try_set!(p_nominal, plant, :V1, 0.0)
    try_set!(p_nominal, plant, :V2, 0.0)
    try_set!(p_nominal, plant, :V13, 0.0)
    try_set!(p_nominal, plant, :V23, 0.0)
    try_set!(p_nominal, plant, :VL1, 0.0)
    try_set!(p_nominal, plant, :VL2, 0.0)
    try_set!(p_nominal, plant, :VL3, 1.0)

    sim = Simulation(
        plant,
        controller,
        copy(h0),
        Dict{Symbol,Any}(
            :t => Float64[],
            :h => Vector{Float64}[],
            :u_Q => Float64[],
            :u_V => Float64[]
        ),
        p_static,
        p_nominal
    )
    return sim
end

function run!(sim::Simulation; duration::Float64, Ts::Float64=1.0)
    steps = Int(duration / Ts)
    t_current = 0.0

    if isempty(sim.history[:t])
        push!(sim.history[:t], t_current)
        push!(sim.history[:h], copy(sim.current_state))
    else
        t_current = sim.history[:t][end]
    end

    for k in 1:steps
        # 1. Compute control
        u_opt = compute_control(sim.controller, sim.current_state, t_current)

        # 2. Record inputs
        push!(sim.history[:u_Q], u_opt.Q1)
        push!(sim.history[:u_V], Float64(u_opt.V1))

        # 3. Update ODE parameters
        p_current_inputs = Dict{Any,Any}()

        # Explicit property access
        sys = sim.plant
        try
            p_current_inputs[sys.Q1] = u_opt.Q1
        catch
        end
        try
            p_current_inputs[sys.Q2] = u_opt.Q2
        catch
        end
        try
            p_current_inputs[sys.V1] = Float64(u_opt.V1)
        catch
        end
        try
            p_current_inputs[sys.V2] = Float64(u_opt.V2)
        catch
        end
        try
            p_current_inputs[sys.V13] = Float64(u_opt.V13)
        catch
        end
        try
            p_current_inputs[sys.V23] = Float64(u_opt.V23)
        catch
        end
        try
            p_current_inputs[sys.VL1] = Float64(u_opt.VL1)
        catch
        end
        try
            p_current_inputs[sys.VL2] = Float64(u_opt.VL2)
        catch
        end
        try
            p_current_inputs[sys.VL3] = Float64(u_opt.VL3)
        catch
        end

        p_step = merge(sim.p_static, sim.p_nominal, p_current_inputs)

        # 4. Integrate
        # Map initial conditions to symbolic state variables
        # Using explicit property access for states h1, h2, h3
        u0_step = [
            sys.h1 => max(sim.current_state[1], 1e-4),
            sys.h2 => max(sim.current_state[2], 1e-4),
            sys.h3 => max(sim.current_state[3], 1e-4)
        ]

        # prob_step = ODEProblem(sim.plant, u0_step, (t_current, t_current + Ts), p_step)
        prob_step = ODEProblem(sim.plant, merge(Dict(u0_step), p_step), (t_current, t_current + Ts))
        sol = solve(prob_step, Rodas5P(), saveat=Ts, reltol=1e-6, abstol=1e-6)

        # Update state
        sim.current_state = [sol[sys.h1][end], sol[sys.h2][end], sol[sys.h3][end]]
        t_current += Ts

        # Record state
        push!(sim.history[:t], t_current)
        push!(sim.history[:h], copy(sim.current_state))
    end
end
