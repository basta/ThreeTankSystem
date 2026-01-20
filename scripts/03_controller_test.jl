using ThreeTankSystem
using Test
using JuMP

function test_mpc_controller()
    println("Initializing Controller...")
    # Initialize controller
    p = TankParameters()
    ctrl = MPCController(params=p, horizon=5, references=[0.5, 0.4, 0.3])

    # Initial state
    h0 = [0.0, 0.0, 0.0]

    println("Step 1 (Should build model)...")
    @time u1 = compute_control(ctrl, h0, 0.0)
    println("Control output 1: ", u1)

    @test !isnothing(ctrl.optimizer_cache)
    model_ref_1 = ctrl.optimizer_cache.model

    # Step 2
    println("Step 2 (Should use cache)...")
    h_next = [0.01, 0.0, 0.0] # Dummy next state
    @time u2 = compute_control(ctrl, h_next, 1.0)
    println("Control output 2: ", u2)

    model_ref_2 = ctrl.optimizer_cache.model

    @test model_ref_1 === model_ref_2
    println("Model reference is identical: ", model_ref_1 === model_ref_2)

    # Step 3: Dynamically fix Q1
    println("Step 3 (Fix Q1=5e-5)...")
    u3 = compute_control(ctrl, h_next, 2.0; fixed_u=Dict(:Q1 => 5e-5))
    println("Control output 3: ", u3)
    @test isapprox(u3.Q1, 5e-5, atol=1e-6)

    # Check internal state
    q1_var = ctrl.optimizer_cache.vars.Q1[1]
    @test is_fixed(q1_var)

    # Step 4: Unfix Q1 (should revert to being optimized)
    println("Step 4 (Unfix Q1)...")
    u4 = compute_control(ctrl, h_next, 3.0; fixed_u=Dict())
    println("Control output 4: ", u4)

    # This assertion is expected to FAIL before the fix
    @test !is_fixed(q1_var)
    println("Q1 is fixed: ", is_fixed(q1_var))

    println("Test passed!")
end

test_mpc_controller()
