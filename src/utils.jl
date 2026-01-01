using GLMakie

# Function to run the visualizer based on your solution
function visualize_tanks(sol, sys)
    t = sol.t
    framerate = 30
    t_interp = range(t[1], t[end], length=300)

    time_index = Observable(1)

    # Lift observables to get current height at the time_index
    # We use sol(t) to interpolate the solution at specific time points
    cur_h1 = @lift([sol(t_interp[$time_index], idxs=sys.h1)])
    cur_h2 = @lift([sol(t_interp[$time_index], idxs=sys.h2)])
    cur_h3 = @lift([sol(t_interp[$time_index], idxs=sys.h3)])

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title=@lift("Time: $(round(t_interp[$time_index], digits=2)) s"),
        xlabel="Tanks", ylabel="Level (m)",
        limits=(0, 4, 0, 0.7))

    # --- Draw Static Elements (Tanks) ---
    # Tank width and positions
    w = 0.5
    pos_x = [1, 3, 2] # Tank 1, 2, and 3 (middle)

    # Tank Outlines (Rectangles)
    # Tank 1
    poly!(ax, Rect(pos_x[1] - w / 2, 0, w, 0.62), color=:transparent, strokecolor=:black, strokewidth=3)
    # Tank 2
    poly!(ax, Rect(pos_x[2] - w / 2, 0, w, 0.62), color=:transparent, strokecolor=:black, strokewidth=3)
    # Tank 3 (Middle)
    poly!(ax, Rect(pos_x[3] - w / 2, 0, w, 0.62), color=:transparent, strokecolor=:black, strokewidth=3)

    # Upper Pipe Connection (Visual reference line at hv=0.3)
    hlines!(ax, [0.3], color=:gray, linestyle=:dash, label="hv (Upper Pipe)")

    # We trick barplot to put bars at specific X coords with specific heights

    # Tank 1 Water
    barplot!(ax, [pos_x[1]], cur_h1, width=w, color=:blue, gap=0)
    # Tank 2 Water
    barplot!(ax, [pos_x[2]], cur_h2, width=w, color=:blue, gap=0)
    # Tank 3 Water
    barplot!(ax, [pos_x[3]], cur_h3, width=w, color=:blue, gap=0)

    # Display the figure
    display(fig)

    # Animation Loop
    for i in 1:length(t_interp)
        time_index[] = i
        sleep(1 / framerate)
    end
end
