using Revise
using ThreeTankSystem
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using JuMP
using HiGHS
using LinearAlgebra

# Simulation Parameters
Ts = 5.0
T_sim = 300.0

tp = TankParameters()
sys = structural_simplify(build_nonlinear_plant(tp))
h0 = [0.0, 0.04, 0.5]
h_ref = [0.5, 0.0, 0.1]

# Controller
lambda_switch = 0.1
settings = InputSettings()
# Enable control for pumps and some valves
set_active!(settings, :Q1, :V1)

# Fix VL3 to always be open (nominal value)
set_nominal!(settings, :VL3, 1.0)
# Other valves (V13, V23, VL1, VL2) default to 0.0, and are NOT active

ctrl = MPCController(params=tp, horizon=5, references=h_ref, lambda_switch=lambda_switch, settings=settings)

# Simulation
sim = Simulation(plant=sys, controller=ctrl, h0=h0)

println("Starting simulation")
run!(sim; duration=T_sim, Ts=Ts)

# Extract results for plotting
t_history = sim.history[:t]
h_history = reduce(vcat, transpose.(sim.history[:h])) # Convert to matrix
u_history_Q = sim.history[:u_Q]
u_history_V = sim.history[:u_V]

# Plotting (reuse previous logic)
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