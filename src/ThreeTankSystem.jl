module ThreeTankSystem

include("types.jl")
include("plant.jl")
include("utils.jl")
include("MLDModel.jl")
include("controllers.jl")
include("simulation.jl")

export build_nonlinear_plant, visualize_tanks, solve_mpc,
    TankParameters, AbstractTankController, MPCController, compute_control,
    Simulation, run!, InputSettings, set_active!, set_nominal!

end # module ThreeTankSystem
