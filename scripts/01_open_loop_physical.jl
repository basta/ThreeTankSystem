using Revise
using ThreeTankSystem
using ModelingToolkit
using OrdinaryDiffEq
using Plots

sys = build_nonlinear_plant()
# 2. Simplify the System
# MTK needs to "compile" the symbolic equations into something numerical.
# It eliminates redundant equations (index reduction) here.
sys_simplified = structural_simplify(sys)

# 3. Define Simulation Conditions
# You can use the symbols directly (sys.x) or names (:x)
tspan = (0.0, 10.0)

sim_cond = [
    sys_simplified.h1 => 0,
    sys_simplified.h2 => 0,
    sys_simplified.h3 => 0,
    sys_simplified.Q1 => 0.003,
    sys_simplified.Q2 => 0.001,
    sys_simplified.VL1 => 0.,
    sys_simplified.VL2 => 0.,
    sys_simplified.VL3 => 1.,
    sys_simplified.V13 => 1.,
    sys_simplified.V23 => 0.,
    sys_simplified.V1 => 1.,
    sys_simplified.V2 => 0.,
]


# 4. Build and Solve
prob = ODEProblem(sys_simplified, sim_cond, tspan)
sol = solve(prob, Tsit5())

visualize_tanks(sol, sys_simplified)

# 5. Plot Results
plot(sol,
    idxs=[sys.h1, sys.h2, sys.h3],
    layout=(3, 1),
    label=["h1" "h2" "h3"],
    title="Open Loop Step Response",
    lw=2
)