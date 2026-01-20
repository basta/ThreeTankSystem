using JuMP, HiGHS

function add_product_constraint!(model, z, delta, f_x)
    # z = delta * f_x
    # delta == 1 => z == f_x
    # delta == 0 => z == 0

    # We use broadcasting for vector inputs
    @constraint(model, [i = eachindex(z)], delta[i] => {z[i] == f_x[i]})
    @constraint(model, [i = eachindex(z)], !delta[i] => {z[i] == 0})
end

function add_switching_cost!(model, x, x_prev_val)
    N = length(x)
    # Create auxiliary variables for absolute difference: d[k] >= |x[k] - x[k-1]|
    d = @variable(model, [1:N], lower_bound = 0)

    @constraint(model, d[1] >= x[1] - x_prev_val)
    @constraint(model, d[1] >= -(x[1] - x_prev_val))

    @constraint(model, [k = 2:N], d[k] >= x[k] - x[k-1])
    @constraint(model, [k = 2:N], d[k] >= -(x[k] - x[k-1]))

    # Return the sum so it can be added to the objective
    return sum(d)
end

function solve_mpc(h0, href, u_prev=Dict(); N=10, p::TankParameters=TankParameters(), fixed_u=Dict())
    mld_model = Model(HiGHS.Optimizer)
    set_silent(mld_model)

    all_inputs = [:Q1, :Q2, :V1, :V2, :V13, :V23, :VL1, :VL2, :VL3]
    optimizing_inputs = setdiff(all_inputs, keys(fixed_u))
    println("Optimizing parameters: ", optimizing_inputs)

    λ_switch = 0.1         # Penalizes switching valve states for V*

    # Use parameters from p
    A = p.A
    az = p.az
    Sh = p.Sh
    g = p.g
    hv = p.hv
    hmax = p.hmax
    Qimax = p.Qimax
    Ts = p.Ts
    SL1 = p.SL1
    SL2 = p.SL2
    SL3 = p.SL3
    SL13 = p.SL13
    SL23 = p.SL23
    S1 = p.S1
    S2 = p.S2

    Q_max = Qimax
    h_max = hmax

    # The inputs
    @variable(mld_model, 0 <= Q1[1:N] <= Q_max)
    @variable(mld_model, 0 <= Q2[1:N] <= Q_max)

    # The connecting valves
    @variable(mld_model, V13[1:N], Bin)
    @variable(mld_model, V23[1:N], Bin)

    # The upper pipe valves
    @variable(mld_model, V1[1:N], Bin)
    @variable(mld_model, V2[1:N], Bin)

    # The leak/output valves 
    @variable(mld_model, VL1[1:N], Bin)
    @variable(mld_model, VL2[1:N], Bin)
    @variable(mld_model, VL3[1:N], Bin)

    # State variables
    @variable(mld_model, 0 <= h1[1:N+1] <= h_max)
    @variable(mld_model, 0 <= h2[1:N+1] <= h_max)
    @variable(mld_model, 0 <= h3[1:N+1] <= h_max)

    # δ0i = 1 if hi > hv, 0 otherwise
    # level threshold indicator, does water level in i reach the upper valve height hv
    @variable(mld_model, δ01[1:N], Bin)
    @variable(mld_model, δ02[1:N], Bin)
    @variable(mld_model, δ03[1:N], Bin)


    # zi3 ≡ Vi3 (hi - h3): flow through the bottom valves. Vi3 must be open
    @variable(mld_model, -0.7 <= z13[1:N] <= 0.7)
    @variable(mld_model, -0.7 <= z23[1:N] <= 0.7)

    # z0i ≡ max{hv, hi} - hv = δ0i * (hi - hv): the height of the water level in i above the upper valve height hv
    @variable(mld_model, -0.7 <= z01[1:N] <= 0.7)
    @variable(mld_model, -0.7 <= z02[1:N] <= 0.7)
    @variable(mld_model, -0.7 <= z03[1:N] <= 0.7)

    # zi ≡ Vi (z0i - z03): flow to the middle. Vi must be open. z0i or z03 must be nonzero; above the upper valve height
    @variable(mld_model, -0.7 <= z1[1:N] <= 0.7)
    @variable(mld_model, -0.7 <= z2[1:N] <= 0.7)

    @variable(mld_model, 0 <= zL1[1:N] <= 0.7) # VL1 * h1
    @variable(mld_model, 0 <= zL2[1:N] <= 0.7) # VL2 * h2
    @variable(mld_model, 0 <= zL3[1:N] <= 0.7) # VL3 * h3

    if haskey(fixed_u, :Q1)
        fix.(Q1, fixed_u[:Q1]; force=true)
    end
    if haskey(fixed_u, :Q2)
        fix.(Q2, fixed_u[:Q2]; force=true)
    end
    if haskey(fixed_u, :V1)
        fix.(V1, fixed_u[:V1]; force=true)
    end
    if haskey(fixed_u, :V2)
        fix.(V2, fixed_u[:V2]; force=true)
    end
    if haskey(fixed_u, :V13)
        fix.(V13, fixed_u[:V13]; force=true)
    end
    if haskey(fixed_u, :V23)
        fix.(V23, fixed_u[:V23]; force=true)
    end
    if haskey(fixed_u, :VL1)
        fix.(VL1, fixed_u[:VL1]; force=true)
    end
    if haskey(fixed_u, :VL2)
        fix.(VL2, fixed_u[:VL2]; force=true)
    end
    if haskey(fixed_u, :VL3)
        fix.(VL3, fixed_u[:VL3]; force=true)
    end

    eps = 1e-5

    for k in 1:N
        # Tank 1
        # δ01[k] == 1 => h1 >= hv
        @constraint(mld_model, δ01[k] => {h1[k] >= hv})
        # δ01[k] == 0 => h1 <= hv - eps
        @constraint(mld_model, !δ01[k] => {h1[k] <= hv - eps})

        # Tank 2
        @constraint(mld_model, δ02[k] => {h2[k] >= hv})
        @constraint(mld_model, !δ02[k] => {h2[k] <= hv - eps})

        # Tank 3
        @constraint(mld_model, δ03[k] => {h3[k] >= hv})
        @constraint(mld_model, !δ03[k] => {h3[k] <= hv - eps})
    end

    # z_i3 = V_i3 * (h_i - h_3)
    add_product_constraint!(mld_model, z13, V13, h1[1:N] .- h3[1:N])
    add_product_constraint!(mld_model, z23, V23, h2[1:N] .- h3[1:N])

    # z_0i = δ0i * (hi - hv)
    add_product_constraint!(mld_model, z01, δ01, h1[1:N] .- hv)
    add_product_constraint!(mld_model, z02, δ02, h2[1:N] .- hv)
    add_product_constraint!(mld_model, z03, δ03, h3[1:N] .- hv)

    # z_i = V_i * (z0i - z03)
    add_product_constraint!(mld_model, z1, V1, z01 .- z03)
    add_product_constraint!(mld_model, z2, V2, z02 .- z03)

    # z = V * h
    add_product_constraint!(mld_model, zL1, VL1, h1[1:N])
    add_product_constraint!(mld_model, zL2, VL2, h2[1:N])
    add_product_constraint!(mld_model, zL3, VL3, h3[1:N])

    coeff = Ts / A


    k_13 = az * SL13 * sqrt(2 * g / hmax)
    k_23 = az * SL23 * sqrt(2 * g / hmax)
    k_1 = az * S1 * sqrt(2 * g / (hmax - hv))
    k_2 = az * S2 * sqrt(2 * g / (hmax - hv))
    k_L1 = az * SL1 * sqrt(2 * g / hmax)
    k_L2 = az * SL2 * sqrt(2 * g / hmax)
    k_L3 = az * SL3 * sqrt(2 * g / hmax)

    fix(h1[1], h0[1]; force=true)
    fix(h2[1], h0[2]; force=true)
    fix(h3[1], h0[3]; force=true)

    @constraint(mld_model, [k = 1:N],
        h1[k+1] == h1[k] + coeff * (Q1[k] - k_1 * z1[k] - k_13 * z13[k] - k_L1 * zL1[k])
    )

    @constraint(mld_model, [k = 1:N],
        h2[k+1] == h2[k] + coeff * (Q2[k] - k_2 * z2[k] - k_23 * z23[k] - k_L2 * zL2[k])
    )

    @constraint(mld_model, [k = 1:N],
        h3[k+1] == h3[k] + coeff * (k_1 * z1[k] + k_2 * z2[k] + k_13 * z13[k] + k_23 * z23[k] - k_L3 * zL3[k])
    )


    @variable(mld_model, 0 <= abs_err_h1[1:N] <= 0.7)
    @variable(mld_model, 0 <= abs_err_h2[1:N] <= 0.7)
    @variable(mld_model, 0 <= abs_err_h3[1:N] <= 0.7)


    switching_penalty_expr = AffExpr(0.0)

    valve_names = [:V1, :V2, :V13, :V23, :VL1, :VL2, :VL3]

    for name in valve_names
        var_ref = mld_model[name]

        val_prev = get(u_prev, name, 0)

        add_to_expression!(switching_penalty_expr,
            add_switching_cost!(mld_model, var_ref, val_prev)
        )
    end



    # Add constraints to force abs_err >= |h - href|
    #    This creates the "V" shape of the absolute value function
    for k in 1:N
        # For h1
        @constraint(mld_model, abs_err_h1[k] >= h1[k+1] - href[1])
        @constraint(mld_model, abs_err_h1[k] >= -(h1[k+1] - href[1]))

        # For h2
        @constraint(mld_model, abs_err_h2[k] >= h2[k+1] - href[2])
        @constraint(mld_model, abs_err_h2[k] >= -(h2[k+1] - href[2]))

        # For h3
        @constraint(mld_model, abs_err_h3[k] >= h3[k+1] - href[3])
        @constraint(mld_model, abs_err_h3[k] >= -(h3[k+1] - href[3]))
    end

    # Note: Q1 and Q2 are already constrained >= 0, so sum(Q) is effectively sum(|Q|)
    @objective(mld_model, Min,
        sum(abs_err_h1[k] + abs_err_h2[k] + abs_err_h3[k] for k in 1:N) +
        1.0 * sum(Q1[k] + Q2[k] for k in 1:N) +
        λ_switch * switching_penalty_expr
    )

    optimize!(mld_model)

    if termination_status(mld_model) != MOI.OPTIMAL
        @warn "MPC Solver Failed or Suboptimal: $(termination_status(mld_model))"
        # Return a safe fallback (everything off/closed)
        return (Q1=0.0, Q2=0.0, V1=0, V2=0, V13=0, V23=0, VL1=0, VL2=0, VL3=0)
    end

    return (
        Q1=value(Q1[1]), Q2=value(Q2[1]),
        V1=round(Int, value(V1[1])), V2=round(Int, value(V2[1])),
        V13=round(Int, value(V13[1])), V23=round(Int, value(V23[1])),
        VL1=round(Int, value(VL1[1])), VL2=round(Int, value(VL2[1])),
        VL3=round(Int, value(VL3[1]))
    )
end