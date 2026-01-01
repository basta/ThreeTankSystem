using Revise
using ThreeTankSystem
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using JuMP
using HiGHS
using LinearAlgebra

# Simulation Parameters
Ts = 1
T_sim = 300.0
steps = Int(T_sim / Ts)

h_ref = [0.5, 0.0, 0.1]

sys = structural_simplify(build_nonlinear_plant())

p_static = Dict(
    sys.A => 0.0154,
    sys.az => 1.0,
    sys.g => 9.81,
    sys.hv => 0.3,
    sys.SL1 => 2e-5, sys.SL2 => 2e-5, sys.SL3 => 2e-5,
    sys.SL13 => 2e-5, sys.SL23 => 2e-5,
    sys.S1 => 2e-5, sys.S2 => 2e-5
)

p_nominal = Dict([
    sys.Q1 => 0,
    sys.Q2 => 0.0,
    sys.V1 => 0.0,
    sys.V2 => 0.0,
    sys.V13 => 0.0,
    sys.V23 => 0.0,
    sys.VL1 => 0.0,
    sys.VL2 => 0.0,
    sys.VL3 => 1.0, # outflow open
])

# :Q1 and :V1 not included, will be optimized
nominal_fixed_u = Dict(
    :Q2 => 0.,
    :V2 => 0.0,
    :V13 => 0.0, :V23 => 0.0,
    :VL1 => 0.0, :VL2 => 0.0,
    :VL3 => 1.0,
)

# initial state
current_h = [0.0, 0.04, 0.5]
u0_plant = [sys.h1 => current_h[1], sys.h2 => current_h[2], sys.h3 => current_h[3]]

h_history = zeros(steps + 1, 3)
u_history_Q = zeros(steps)
u_history_V = zeros(steps)
t_history = 0:Ts:T_sim

prob = ODEProblem(sys, u0_plant, (0.0, T_sim), p_nominal)

println("Starting simulation")
last_u = Dict()
for k in 1:steps
    global current_h, h_ref, last_u

    u_opt = solve_mpc(current_h, h_ref, last_u; fixed_u=nominal_fixed_u, N=5, Ts=Ts)
    last_u = u_opt
    u_history_Q[k] = u_opt.Q1
    u_history_V[k] = u_opt.V1
    h_history[k, :] = current_h

    p_current_inputs = Dict(
        sys.Q1 => u_opt.Q1,
        sys.V1 => Float64(u_opt.V1),
        sys.V2 => Float64(u_opt.V2),
        sys.Q2 => u_opt.Q2,
        # Ensure other inputs are maintained from nominal if not optimized
        sys.VL3 => p_nominal[sys.VL3],
        sys.V13 => p_nominal[sys.V13],
        sys.V23 => p_nominal[sys.V23],
        sys.VL1 => p_nominal[sys.VL1],
        sys.VL2 => p_nominal[sys.VL2]
    )

    p_step = merge(p_static, p_current_inputs)

    u0_step = [
        sys.h1 => max(current_h[1], 1e-4),
        sys.h2 => max(current_h[2], 1e-4),
        sys.h3 => max(current_h[3], 1e-4)
    ]

    prob_step = ODEProblem(sys, u0_step, ((k - 1) * Ts, k * Ts), p_step)

    sol = solve(prob_step, Rodas5P(), saveat=Ts, reltol=1e-6, abstol=1e-6)
    # Explicitly extract states in the correct order [h1, h2, h3]
    current_h = [sol[sys.h1][end], sol[sys.h2][end], sol[sys.h3][end]]
end
h_history[end, :] = current_h

p1 = plot(t_history, h_history[:, 1], label="h1", lw=2, color=:blue,
    ylabel="Level [m]", title="Nominal Operation", legend=:bottomright)
plot!(p1, t_history, h_history[:, 2], label="h2", lw=2, color=:green)
plot!(p1, t_history, h_history[:, 3], label="h3", lw=2, color=:red)
hline!(p1, [h_ref[1]], label="Ref h1", linestyle=:dash, color=:blue)
hline!(p1, [h_ref[2]], label="Ref h2", linestyle=:dash, color=:green)
hline!(p1, [h_ref[3]], label="Ref h3", linestyle=:dash, color=:red)

p2 = plot(t_history[1:end-1], u_history_Q, label="Q1 (Pump)", linetype=:steppost,
    ylabel="Flow", color=:black)
p3 = plot(t_history[1:end-1], u_history_V, label="V1 (Valve)", linetype=:steppost,
    ylabel="Binary", xlabel="Time [s]", color=:green)

plot(p1, p2, p3, layout=(3, 1), size=(600, 800))