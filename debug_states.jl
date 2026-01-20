using ThreeTankSystem
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Setup
tp = TankParameters()
sys = structural_simplify(build_nonlinear_plant(tp))

println("Unknowns: ", unknowns(sys))

# Check properties
println("sys.h1: ", sys.h1)
println("sys.h2: ", sys.h2)
println("sys.h3: ", sys.h3)

# Check equality
println("sys.h1 == unknowns(sys)[3]? ", isequal(sys.h1, unknowns(sys)[3]))
println("sys.h2 == unknowns(sys)[2]? ", isequal(sys.h2, unknowns(sys)[2]))
println("sys.h3 == unknowns(sys)[1]? ", isequal(sys.h3, unknowns(sys)[1]))

# Simulate one step to check order
using OrdinaryDiffEq
u0 = [sys.h1 => 1.0, sys.h2 => 2.0, sys.h3 => 3.0]
prob = ODEProblem(sys, u0, (0.0, 0.1), [])
sol = solve(prob, Tsit5())

println("Sol at end: ", sol[end])
println("Sol[sys.h1]: ", sol[sys.h1][end])
println("Sol[sys.h2]: ", sol[sys.h2][end])
println("Sol[sys.h3]: ", sol[sys.h3][end])

extracted = [sol[sys.h1][end], sol[sys.h2][end], sol[sys.h3][end]]
println("Extracted: ", extracted)
