### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 51116c5a-f643-11f0-a92d-2d83747e0a66
begin
    using Pkg
    Pkg.activate("..") # Activate the local ThreeTankSystem environment
    using PlutoUI
    using Plots
    using ModelingToolkit
    using OrdinaryDiffEq
    using JuMP
    using HiGHS
    using LinearAlgebra
    using Printf
    
    # Load the local package
    using ThreeTankSystem
    
    # Set plotting theme
    theme(:wong)
end

# ╔═╡ fd4b345b-1191-47bd-8d6b-39bbc88143fc
function draw_system_schematic()
    # Helper to draw a tank
    tank(x, y, label) = Shape([x, x+1.2, x+1.2, x], [y, y, y+3, y+3])
    
    # Helper to draw a valve (bowtie shape)
    function valve(x, y, label_text, label_offset=:top)
        w, h = 0.25, 0.15
        # Bowtie polygon
        v_shape = Shape([x-w, x+w, x+w, x-w], [y+h, y-h, y+h, y-h])
        
        # Determine label position
        ly = label_offset == :top ? y + 0.3 : y - 0.4
        return (v_shape, (x, ly, text(label_text, 9, :black)))
    end

    p = plot(aspect_ratio=:equal, axis=nothing, border=:none, size=(800, 400))
    
    # --- TANKS ---
    # T1 (Left), T3 (Middle), T2 (Right) - Coordinates shifted for spacing
    plot!(p, tank(0, 0, ""), color=:white, linecolor=:black, lw=2, label="")
    annotate!(0.6, 1.5, text("1", 20, :black))
    
    plot!(p, tank(3, 0, ""), color=:white, linecolor=:black, lw=2, label="")
    annotate!(3.6, 1.5, text("3", 20, :black))
    
    plot!(p, tank(6, 0, ""), color=:white, linecolor=:black, lw=2, label="")
    annotate!(6.6, 1.5, text("2", 20, :black))

    # --- PUMPS ---
    # Pump 1 (feeds T1)
    plot!(p, [-0.5, 0.6, 0.6], [3.5, 3.5, 3.0], color=:black, lw=2, label="") # Pipe
    scatter!(p, [-0.5], [3.5], markersize=12, color=:white, markerstrokecolor=:black, label="") # Pump symbol
    annotate!(-0.5, 4.0, text("Q1", 10, :bold))
    
    # Pump 2 (feeds T2)
    plot!(p, [7.7, 6.6, 6.6], [3.5, 3.5, 3.0], color=:black, lw=2, label="") # Pipe
    scatter!(p, [7.7], [3.5], markersize=12, color=:white, markerstrokecolor=:black, label="")
    annotate!(7.7, 4.0, text("Q2", 10, :bold))

    # --- CONNECTIONS ---
    # Upper Pipes (at h_v ~ 2.0)
    plot!(p, [1.2, 3.0], [2.0, 2.0], color=:gray, lw=4, alpha=0.5, label="") # T1-T3
    plot!(p, [4.2, 6.0], [2.0, 2.0], color=:gray, lw=4, alpha=0.5, label="") # T3-T2
    
    # Lower Pipes (at bottom ~ 0.5)
    plot!(p, [1.2, 3.0], [0.5, 0.5], color=:gray, lw=4, alpha=0.5, label="") # T1-T3
    plot!(p, [4.2, 6.0], [0.5, 0.5], color=:gray, lw=4, alpha=0.5, label="") # T3-T2

    # --- VALVES ---
    # Draw valves and collect annotations
    # Upper V1, V2
    v1, t1 = valve(2.1, 2.0, "V1")
    v2, t2 = valve(5.1, 2.0, "V2")
    
    # Lower V13, V23
    v13, t13 = valve(2.1, 0.5, "V13", :bottom)
    v23, t23 = valve(5.1, 0.5, "V23", :bottom)
    
    # Drains VL1, VL3, VL2
    plot!(p, [0.6, 0.6], [0, -0.8], color=:black, lw=2, label="") # Pipe T1
    vl1, tl1 = valve(0.6, -0.5, "VL1", :left)
    annotate!(0.9, -0.5, "QL1")

    plot!(p, [3.6, 3.6], [0, -0.8], color=:black, lw=2, label="") # Pipe T3
    vl3, tl3 = valve(3.6, -0.5, "VL3", :left)
    annotate!(3.9, -0.5, "QN3")

    plot!(p, [6.6, 6.6], [0, -0.8], color=:black, lw=2, label="") # Pipe T2
    vl2, tl2 = valve(6.6, -0.5, "VL2", :left)
    annotate!(6.9, -0.5, "QL2")

    # Render Valves
    for v in [v1, v2, v13, v23, vl1, vl3, vl2]
        plot!(p, v, color=:white, linecolor=:black, lw=1, label="")
    end
    
    # Render Text Labels
    for (tx, ty, tstr) in [t1, t2, t13, t23]
        annotate!(tx, ty, tstr)
    end
    
    # Height Labels
    plot!(p, [-0.2, -0.2], [0, 2], arrow=:both, color=:black, lw=1, label="")
    annotate!(-0.5, 1.0, "hv")

    return p
end

# ╔═╡ 5ead0ba0-20c0-4f17-b051-993f9718e9e2
md"""
# COSY Three Tank Benchmark System
**Project:** BP33 Cybernetics and Robotics  
**Topic:** Hybrid Model Predictive Control

## 1. System Topology
The system consists of three tanks interacting via pipes at two different levels. 
* **Inputs:** Two pumps $Q_1, Q_2$ feed Tank 1 and Tank 2 respectively.
* **Interconnections:**
    * **Upper ($h > h_v$):** Controlled by valves $V_1$ (T1-T3) and $V_2$ (T2-T3).
    * **Lower ($h > 0$):** Controlled by valves $V_{13}$ (T1-T3) and $V_{23}$ (T2-T3).
* **Drains:** All three tanks have drain valves $V_{L1}, V_{L3}, V_{L2}$ (represented as leaks or controlled outflows).

## 2. Differential Equations
The dynamics are derived from Torricelli's Law. The central Tank 3 acts as a mixing vessel fed by the outer tanks.

$$
\begin{aligned}
A \dot{h}_1 &= Q_1 - \underbrace{Q_{13}^{(V1)}}_{\text{Upper}} - \underbrace{Q_{13}^{(V13)}}_{\text{Lower}} - Q_{L1} \\
A \dot{h}_2 &= Q_2 - \underbrace{Q_{23}^{(V2)}}_{\text{Upper}} - \underbrace{Q_{23}^{(V23)}}_{\text{Lower}} - Q_{L2} \\
A \dot{h}_3 &= \underbrace{(Q_{13}^{(V1)} + Q_{13}^{(V13)})}_{\text{In from T1}} + \underbrace{(Q_{23}^{(V2)} + Q_{23}^{(V23)})}_{\text{In from T2}} - Q_{N3}
\end{aligned}
$$

The flow $Q_{ij}$ through any valve depends on the valve state ($V \in \{0,1\}$) and the square root of the level difference:
$$Q_{ij} = V \cdot a_{z} S \text{sign}(h_i - h_j) \sqrt{2g |h_i - h_j|}$$
"""

# ╔═╡ 24f138d6-c448-4c12-b2bf-d4b0de67fda9
draw_system_schematic()

# ╔═╡ 884615d1-9dbe-4ec4-9394-d4dff839ccd2
function run_simulation_scenario(
    h0_in, ref_in, horizon_in, lambda_in, duration_in,
    # Fault Flags
    drain_blocked, v1_stuck_closed
)
    # 1. Setup Physical Plant
    tp = TankParameters()
    plant_sys = structural_simplify(build_nonlinear_plant(tp))
    
    # 2. Configure Controller & Faults
    settings = InputSettings()
    
    # Standard Setup: Q1 and V1 are active controls
    # V2 is usually controlled too, let's enable it for full control
    set_active!(settings, :Q1, :V1, :V2, :VL3) 
    
    # --- FAULT INJECTION ---
    # If drain is blocked, we tell the controller "VL3 is fixed to 0" 
    # (reconfiguration) OR we force the plant input to 0 (robustness).
    # Here we simulate 'Constraints': The controller knows it cannot use these actuators.
    
    if drain_blocked
        set_fixed!(settings, :VL3, true)
        set_nominal!(settings, :VL3, 0.0)
    end
    
    if v1_stuck_closed
        set_fixed!(settings, :V1, true)
        set_nominal!(settings, :V1, 0.0)
    end

    # 3. Build MPC
    ctrl = MPCController(
        params=tp, 
        horizon=horizon_in, 
        references=ref_in, 
        lambda_switch=lambda_in, 
        settings=settings
    )

    # 4. Run
    sim = Simulation(plant=plant_sys, controller=ctrl, h0=h0_in)
    run!(sim; duration=duration_in, Ts=5.0) # Ts from your types.jl defaults
    
    return sim
end

# ╔═╡ 1f084c9c-fe0a-4295-9472-baf8195c5015
function run_simulation_matrix(
    h0_in, ref_in, horizon_in, lambda_in, duration_in,
    # Control Flags (Booleans)
    use_Q1, use_Q2,
    use_V1, use_V2,
    use_V13, use_V23,
    use_VL1, use_VL2, use_VL3,
    # Special Nominal Override
    vl3_normally_open
)
    # 1. Setup Physical Plant
    tp = TankParameters()
    plant_sys = structural_simplify(build_nonlinear_plant(tp))
    
    # 2. Configure Controller Inputs
    settings = InputSettings()
    
    # We build the set of active inputs dynamically
    if use_Q1 push!(settings.enabled_inputs, :Q1) end
    if use_Q2 push!(settings.enabled_inputs, :Q2) end
    
    if use_V1 push!(settings.enabled_inputs, :V1) end
    if use_V2 push!(settings.enabled_inputs, :V2) end
    
    if use_V13 push!(settings.enabled_inputs, :V13) end
    if use_V23 push!(settings.enabled_inputs, :V23) end
    
    if use_VL1 push!(settings.enabled_inputs, :VL1) end
    if use_VL2 push!(settings.enabled_inputs, :VL2) end
    if use_VL3 push!(settings.enabled_inputs, :VL3) end

    # 3. Handle Nominal Values (When NOT optimized)
    # By default in InputSettings, everything is 0.0 (Closed/Off)
    # We allow one special override for VL3 because it is often a passive drain.
    if !use_VL3 && vl3_normally_open
        set_nominal!(settings, :VL3, 1.0)
    else
        # Ensure it is 0 if not used (redundant but safe)
        set_nominal!(settings, :VL3, 0.0)
    end


    # 4. Build MPC
    ctrl = MPCController(
        params=tp, 
        horizon=horizon_in, 
        references=ref_in, 
        lambda_switch=lambda_in, 
        settings=settings
    )

    # 5. Run
    sim = Simulation(plant=plant_sys, controller=ctrl, h0=h0_in)
    run!(sim; duration=duration_in, Ts=5.0)
    
    return sim
end

# ╔═╡ 54d037a4-3bff-429f-8e75-5a9334545602
begin
    md"""
    ### Simulation Control
    
    **Initial Levels ($h_0$):** T1: $(@bind h1_0 Slider(0.0:0.05:0.6; default=0.0, show_value=true)) m  
    T3: $(@bind h3_0 Slider(0.0:0.05:0.6; default=0.2, show_value=true)) m  
    T2: $(@bind h2_0 Slider(0.0:0.05:0.6; default=0.0, show_value=true)) m
    
    **Target Levels ($h_{ref}$):** Ref 1: $(@bind ref1 Slider(0.0:0.05:0.6; default=0.4, show_value=true))  
    Ref 3: $(@bind ref3 Slider(0.0:0.05:0.6; default=0.2, show_value=true))  
    Ref 2: $(@bind ref2 Slider(0.0:0.05:0.6; default=0.1, show_value=true))

    ---
    Stability/Switching ($\lambda$): $(@bind lambda_mpc Slider(0.0:0.01:1.0; default=0.1, show_value=true))
    
    ---
    ### Actuator Selection
    *Select which controls the MPC is allowed to use. Unchecked items remain **Closed (0)**.*
    
    **Pumps:** $(@bind ctrl_Q1 CheckBox(default=true)) **Q1** (Pump 1)  
    $(@bind ctrl_Q2 CheckBox(default=true)) **Q2** (Pump 2)
    
    **Connecting Valves:** $(@bind ctrl_V1 CheckBox(default=true)) **V1** (Upper 1-3)  
    $(@bind ctrl_V2 CheckBox(default=true)) **V2** (Upper 2-3)  
    $(@bind ctrl_V13 CheckBox(default=false)) **V13** (Lower 1-3)  
    $(@bind ctrl_V23 CheckBox(default=false)) **V23** (Lower 2-3)
    
    **Drains:** $(@bind ctrl_VL1 CheckBox(default=false)) **VL1** (Drain 1)  
    $(@bind ctrl_VL2 CheckBox(default=false)) **VL2** (Drain 2)  
    $(@bind ctrl_VL3 CheckBox(default=false)) **VL3** (Drain 3)
    
    *Special:* $(@bind nom_VL3 CheckBox(default=true)) **Keep VL3 Open** if not controlled?
    
    ---
    $(@bind run_btn Button("▶️ RUN SCENARIO"))
    """
end

# ╔═╡ 77deb1da-cd66-4e28-8749-3914865bb0ff
begin
    run_btn # Dependency on button
	N_mpc = 5
    
    # Prepare Inputs
    local_h0 = [h1_0, h2_0, h3_0] 
    local_ref = [ref1, ref2, ref3]
    
    # Execute
    sim_result = run_simulation_matrix(
        local_h0, local_ref, N_mpc, lambda_mpc, 300.0,
        ctrl_Q1, ctrl_Q2,
        ctrl_V1, ctrl_V2,
        ctrl_V13, ctrl_V23,
        ctrl_VL1, ctrl_VL2, ctrl_VL3,
        nom_VL3
    );
    
end

# ╔═╡ d5c9e170-e0f9-465f-88af-70976b17c093
begin
    # Data Extraction
    t = sim_result.history[:t]
    h_data = reduce(vcat, transpose.(sim_result.history[:h])) # Matrix [N x 3]
    u_Q = sim_result.history[:u_Q]
    u_V = sim_result.history[:u_V] # This is likely a vector of vectors or matrix
    
    # 1. Main Plot: Water Levels
    p_levels = plot(t, h_data, label=["T1" "T2" "T3"], lw=3,
        xlabel="Time [s]", ylabel="Level [m]", title="Water Level Response",
        color=[:blue :green :orange], legend=:topright)
        
    # Add Setpoints
    hline!(p_levels, [ref1], style=:dash, color=:blue, label="Ref 1", alpha=0.5)
    hline!(p_levels, [ref2], style=:dash, color=:green, label="Ref 2", alpha=0.5)
    hline!(p_levels, [ref3], style=:dash, color=:orange, label="Ref 3", alpha=0.5)
    
    # 2. Actuator Plot (Pump)
    p_pump = plot(t[1:end-1], u_Q, linetype=:steppost, fill=(0, 0.2, :blue),
        ylabel="Flow [m3/s]", title="Pump Q1 Input", legend=false, color=:blue)
        
    # 3. Valve States (Heatmap or Step)
    # u_V is likely [V1, V2, VL3] or similar. Let's process it.
    # We transpose to get [Time x Valves]
    V_mat = reduce(vcat, transpose.(u_V))
    p_valves = plot(t[1:end-1], V_mat, linetype=:steppost,
        label=["V1" "V2" "VL3"], title="Discrete Valve States",
        yticks=[0, 1], ylabel="Open/Closed", lw=2)

    # Combine
    plot(p_levels, p_pump, p_valves, layout=@layout([a; b c]), size=(800, 800))
end

# ╔═╡ a139df4c-958e-41e8-96fd-450a312843f5


# ╔═╡ 62adcbe2-5990-489f-aadd-791ff9846288
function draw_dynamic_frame(h_current, u_current_Q, u_current_V)
    # Re-use the static background logic, but we will overlay dynamic elements
    # We copy the setup from your draw_system_schematic for consistency
    
    p = plot(aspect_ratio=:equal, axis=nothing, border=:none, size=(800, 400))
    
    # --- CONSTANTS & SCALING ---
    # Map physical height (e.g., 0.6m) to visual height (3.0 units)
    H_SCALE = 3.0 / 0.6
    
    # Valve State Colors (Green = Open, White = Closed)
    v_color(is_open) = is_open > 0.5 ? :green : :white
    
    # --- 1. DRAW TANKS (Outlines) ---
    tank_xs = [0.0, 6.0, 3.0] # T1, T2, T3 (Note T3 is middle visually, but index 3)
    tank_labels = ["1", "2", "3"]
    
    for (i, x) in enumerate(tank_xs)
        # Draw Tank Outline
        plot!(p, Shape([x, x+1.2, x+1.2, x], [0, 0, 3, 3]), 
              color=:white, linecolor=:black, lw=2, label="")
        
        # --- DYNAMIC WATER ---
        # Draw Water Level Rect (Blue)
        # Clamp height to visual max (3.0) to prevent drawing spillover
        h_visual = clamp(h_current[i] * H_SCALE, 0, 3.0)
        if h_visual > 0.01
            plot!(p, Shape([x, x+1.2, x+1.2, x], [0, 0, h_visual, h_visual]), 
                  color=RGBA(0, 0, 1, 0.5), linecolor=nothing, label="")
        end
        
        annotate!(x+0.6, 1.5, text(tank_labels[i], 20, :black))
    end

    # --- 2. DYNAMIC PUMPS ---
    # Pump 1 (Left)
    q1_active = u_current_Q[1] > 0
    c_p1 = q1_active ? :green : :white
    scatter!(p, [-0.5], [3.5], markersize=12, color=c_p1, markerstrokecolor=:black, label="")
    annotate!(-0.5, 4.0, text("Q1: $(round(u_current_Q[1], digits=3))", 8, :black))

    # Pump 2 (Right)
    # Check if Q2 exists in input vector (it might be index 2)
    q2_val = length(u_current_Q) >= 2 ? u_current_Q[2] : 0.0
    c_p2 = q2_val > 0 ? :green : :white
    scatter!(p, [7.7], [3.5], markersize=12, color=c_p2, markerstrokecolor=:black, label="")
    annotate!(7.7, 4.0, text("Q2: $(round(q2_val, digits=3))", 8, :black))

    # --- 3. DYNAMIC VALVES ---
    # Helper to draw valve bowtie with dynamic color
    function draw_valve(x, y, label_text, is_open)
        w, h = 0.25, 0.15
        v_shape = Shape([x-w, x+w, x+w, x-w], [y+h, y-h, y+h, y-h])
        plot!(p, v_shape, color=v_color(is_open), linecolor=:black, lw=1, label="")
        annotate!(x, y+0.3, text(label_text, 8, :black))
    end
    
    # Mapping simulation inputs to valve positions (Logic depends on how MPC output is ordered)
    # Assuming u_current_V order: [V1, V2, V13, V23, VL1, VL2, VL3] based on typical setup
    # If using your checkbox logic, we need to be careful. 
    # For visualization, let's map by name if possible, or assume a standard full vector.
    # Here we assume the full vector is passed or we default to 0.
    
    # Extract values safely (default to 0 if index out of bounds)
    get_v(idx) = length(u_current_V) >= idx ? u_current_V[idx] : 0.0
    
    # Upper
    draw_valve(2.1, 2.0, "V1", get_v(1))
    draw_valve(5.1, 2.0, "V2", get_v(2))
    
    # Lower (if they exist in vector, usually indices 3,4)
    draw_valve(2.1, 0.5, "V13", get_v(3))
    draw_valve(5.1, 0.5, "V23", get_v(4))
    
    # Drains (indices 5,6,7)
    draw_valve(0.6, -0.5, "VL1", get_v(5)) # T1 Drain
    draw_valve(6.6, -0.5, "VL2", get_v(6)) # T2 Drain
    draw_valve(3.6, -0.5, "VL3", get_v(7)) # T3 Drain

    # Pipes (Static gray lines for context)
    plot!(p, [1.2, 3.0], [2.0, 2.0], color=:gray, lw=4, alpha=0.3, label="")
    plot!(p, [4.2, 6.0], [2.0, 2.0], color=:gray, lw=4, alpha=0.3, label="")
    plot!(p, [1.2, 3.0], [0.5, 0.5], color=:gray, lw=4, alpha=0.3, label="")
    plot!(p, [4.2, 6.0], [0.5, 0.5], color=:gray, lw=4, alpha=0.3, label="")

    return p
end

# ╔═╡ 115647e0-567a-45e6-ad6f-64284aee59af
@bind time_idx Slider(1:1:(length(sim_result.history[:t])-1); show_value=false)

# ╔═╡ 1d19d08b-b7a2-4c5f-bd66-31d8d5a936bc
begin
    # 1. Setup Data
    # Extract time and state matrices
    anim_t = sim_result.history[:t]
    anim_h = reduce(vcat, transpose.(sim_result.history[:h])) # [N x 3]
    
    # u_Q and u_V are likely Vector of Vectors in DifferentialEquations solution style
    # We convert them to matrices for easier indexing
    anim_uQ = sim_result.history[:u_Q] 
    anim_uV = sim_result.history[:u_V]
    
    max_idx = length(anim_t)
    
    md"""
    
    Current Time: **$(round(anim_t[time_idx], digits=1)) s**
    
    $(
        # Call the drawing function with data at current index
        begin
            # Extract snapshots
            h_snap = anim_h[time_idx, :]
            Q_snap = anim_uQ[time_idx] # Vector
            V_snap = anim_uV[time_idx] # Vector
            
            draw_dynamic_frame(h_snap, Q_snap, V_snap)
        end
    )
    """
end

# ╔═╡ fab83415-032b-4a42-80cc-84d2b41a546f


# ╔═╡ Cell order:
# ╠═51116c5a-f643-11f0-a92d-2d83747e0a66
# ╟─fd4b345b-1191-47bd-8d6b-39bbc88143fc
# ╟─5ead0ba0-20c0-4f17-b051-993f9718e9e2
# ╠═24f138d6-c448-4c12-b2bf-d4b0de67fda9
# ╟─884615d1-9dbe-4ec4-9394-d4dff839ccd2
# ╠═1f084c9c-fe0a-4295-9472-baf8195c5015
# ╟─54d037a4-3bff-429f-8e75-5a9334545602
# ╠═77deb1da-cd66-4e28-8749-3914865bb0ff
# ╟─d5c9e170-e0f9-465f-88af-70976b17c093
# ╠═a139df4c-958e-41e8-96fd-450a312843f5
# ╟─62adcbe2-5990-489f-aadd-791ff9846288
# ╠═115647e0-567a-45e6-ad6f-64284aee59af
# ╠═1d19d08b-b7a2-4c5f-bd66-31d8d5a936bc
# ╠═fab83415-032b-4a42-80cc-84d2b41a546f
