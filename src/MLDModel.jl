using JuMP, HiGHS

struct MPCModelRefs
    model::Model
    # Variables that constitute the inputs (for fixing/extracting)
    vars::NamedTuple
    # References for updating parameters
    h0_vars::Vector{VariableRef}
    href_constraints::Vector{Vector{ConstraintRef}} # [tank_idx][step] => 2 constraints
    switch_constraints::Dict{Symbol,Vector{ConstraintRef}} # input_name => 2 constraints (for step 1)
end

function add_product_constraint!(model, z, delta, f_x)
    @constraint(model, [i = eachindex(z)], delta[i] => {z[i] == f_x[i]})
    @constraint(model, [i = eachindex(z)], !delta[i] => {z[i] == 0})
end

function add_switching_cost!(model, x, x_prev_val, refs_store)
    N = length(x)
    # Create auxiliary variables for absolute difference: d[k] >= |x[k] - x[k-1]|
    d = @variable(model, [1:N], lower_bound = 0)

    # d[1] >= x[1] - x_prev_val  => d[1] - x[1] >= -x_prev_val
    c1 = @constraint(model, d[1] - x[1] >= -x_prev_val)
    # d[1] >= -(x[1] - x_prev_val) => d[1] + x[1] >= x_prev_val
    c2 = @constraint(model, d[1] + x[1] >= x_prev_val)

    # Store these to update x_prev_val later
    push!(refs_store, c1)
    push!(refs_store, c2)

    @constraint(model, [k = 2:N], d[k] >= x[k] - x[k-1])
    @constraint(model, [k = 2:N], d[k] >= -(x[k] - x[k-1]))

    return sum(d)
end

function build_mpc_model(p::TankParameters, N::Int, λ_switch::Float64)
    mld_model = Model(HiGHS.Optimizer)
    set_silent(mld_model)



    # Unpack parameters
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

    # --- Variables ---
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

    @constraint(mld_model, [k = 1:N],
        h1[k+1] == h1[k] + coeff * (Q1[k] - k_1 * z1[k] - k_13 * z13[k] - k_L1 * zL1[k])
    )

    @constraint(mld_model, [k = 1:N],
        h2[k+1] == h2[k] + coeff * (Q2[k] - k_2 * z2[k] - k_23 * z23[k] - k_L2 * zL2[k])
    )

    @constraint(mld_model, [k = 1:N],
        h3[k+1] == h3[k] + coeff * (k_1 * z1[k] + k_2 * z2[k] + k_13 * z13[k] + k_23 * z23[k] - k_L3 * zL3[k])
    )

    # --- Objective Terms ---
    @variable(mld_model, 0 <= abs_err_h1[1:N] <= 0.7)
    @variable(mld_model, 0 <= abs_err_h2[1:N] <= 0.7)
    @variable(mld_model, 0 <= abs_err_h3[1:N] <= 0.7)


    switching_penalty_expr = AffExpr(0.0)

    switch_constraints = Dict{Symbol,Vector{ConstraintRef}}()
    vars = (
        Q1=Q1, Q2=Q2, V1=V1, V2=V2, V13=V13, V23=V23, VL1=VL1, VL2=VL2, VL3=VL3,
        h1=h1, h2=h2, h3=h3
    )

    valve_names = [:V1, :V2, :V13, :V23, :VL1, :VL2, :VL3]

    for name in valve_names
        var_ref = vars[name]
        sw_refs = Vector{ConstraintRef}()
        cost = add_switching_cost!(mld_model, var_ref, 0.0, sw_refs) # Init with 0.0
        add_to_expression!(switching_penalty_expr, cost)
        switch_constraints[name] = sw_refs
    end

    # Reference tracking constraints
    # abs_err >= h - href => abs_err - h >= -href
    # abs_err >= -(h - href) => abs_err + h >= href
    href_constraints = [Vector{ConstraintRef}(undef, 2 * N) for _ in 1:3]

    for k in 1:N
        # h1
        c1 = @constraint(mld_model, abs_err_h1[k] - h1[k+1] >= -0.0) # Init href=0
        c2 = @constraint(mld_model, abs_err_h1[k] + h1[k+1] >= 0.0)
        href_constraints[1][2*k-1] = c1
        href_constraints[1][2*k] = c2

        # h2
        c3 = @constraint(mld_model, abs_err_h2[k] - h2[k+1] >= -0.0)
        c4 = @constraint(mld_model, abs_err_h2[k] + h2[k+1] >= 0.0)
        href_constraints[2][2*k-1] = c3
        href_constraints[2][2*k] = c4

        # h3
        c5 = @constraint(mld_model, abs_err_h3[k] - h3[k+1] >= -0.0)
        c6 = @constraint(mld_model, abs_err_h3[k] + h3[k+1] >= 0.0)
        href_constraints[3][2*k-1] = c5
        href_constraints[3][2*k] = c6
    end

    @objective(mld_model, Min,
        sum(abs_err_h1[k] + abs_err_h2[k] + abs_err_h3[k] for k in 1:N) +
        1.0 * sum(Q1[k] + Q2[k] for k in 1:N) +
        λ_switch * switching_penalty_expr
    )

    mpc_refs = MPCModelRefs(
        mld_model,
        vars,
        [h1[1], h2[1], h3[1]],
        href_constraints,
        switch_constraints
    )

    return mpc_refs
end

function update_mpc_model!(mpc_refs::MPCModelRefs, h0, href, u_prev=Dict(), fixed_u=Dict())
    model = mpc_refs.model

    # Update initial state
    fix(mpc_refs.h0_vars[1], h0[1]; force=true)
    fix(mpc_refs.h0_vars[2], h0[2]; force=true)
    fix(mpc_refs.h0_vars[3], h0[3]; force=true)

    # Update reference constraints
    N = length(mpc_refs.href_constraints[1]) ÷ 2
    for tank_idx in 1:3
        ref_val = href[tank_idx]
        constrs = mpc_refs.href_constraints[tank_idx]
        for k in 1:N
            c1 = constrs[2*k-1] # abs_err - h >= -href
            c2 = constrs[2*k]   # abs_err + h >= href
            set_normalized_rhs(c1, -ref_val)
            set_normalized_rhs(c2, ref_val)
        end
    end

    # Update switching cost constraints
    # d[1] - x[1] >= -x_prev_val
    # d[1] + x[1] >= x_prev_val
    for (name, constrs) in mpc_refs.switch_constraints
        val_prev = get(u_prev, name, 0.0)
        # Check if val_prev is boolean/int, convert to float
        val_float = Float64(val_prev)
        set_normalized_rhs(constrs[1], -val_float)
        set_normalized_rhs(constrs[2], val_float)
    end

    #  Fix inputs
    for (name, val) in fixed_u
        if haskey(mpc_refs.vars, name)
            fix.(mpc_refs.vars[name], val; force=true)
        end
    end
end

function extract_solution(mpc_refs::MPCModelRefs)
    model = mpc_refs.model
    if termination_status(model) != MOI.OPTIMAL
        @warn "MPC Solver Failed or Suboptimal: $(termination_status(model))"
        # Return fallback
        return (Q1=0.0, Q2=0.0, V1=0, V2=0, V13=0, V23=0, VL1=0, VL2=0, VL3=0)
    end

    v = mpc_refs.vars
    return (
        Q1=value(v.Q1[1]), Q2=value(v.Q2[1]),
        V1=round(Int, value(v.V1[1])), V2=round(Int, value(v.V2[1])),
        V13=round(Int, value(v.V13[1])), V23=round(Int, value(v.V23[1])),
        VL1=round(Int, value(v.VL1[1])), VL2=round(Int, value(v.VL2[1])),
        VL3=round(Int, value(v.VL3[1]))
    )
end

function solve_mpc(h0, href, u_prev=Dict(); N=10, p::TankParameters=TankParameters(), fixed_u=Dict(), λ_switch=0.1)
    # Convenience function using a throw-away model
    mpc_refs = build_mpc_model(p, N, λ_switch)
    update_mpc_model!(mpc_refs, h0, href, u_prev, fixed_u)
    optimize!(mpc_refs.model)
    return extract_solution(mpc_refs)
end