Base.@kwdef struct TankParameters
    A::Float64 = 0.0154
    az::Float64 = 1.0
    g::Float64 = 9.81
    Sh::Float64 = 2e-5
    hv::Float64 = 0.3
    hmax::Float64 = 0.62
    Qimax::Float64 = 1e-4
    Ts::Float64 = 5.0
    SL1::Float64 = 2e-5
    SL2::Float64 = 2e-5
    SL3::Float64 = 2e-5
    SL13::Float64 = 2e-5
    SL23::Float64 = 2e-5
    S1::Float64 = 2e-5
    S2::Float64 = 2e-5
end

Base.@kwdef mutable struct InputSettings
    # Set of symbols that the MPC can control
    # e.g. [:Q1, :V1, :VL3]
    enabled_inputs::Set{Symbol} = Set{Symbol}()

    # Default values for inputs when they are NOT enabled
    nominal_values::Dict{Symbol,Float64} = Dict(
        :Q1 => 0.0, :Q2 => 0.0,
        :V1 => 0.0, :V2 => 0.0,
        :V13 => 0.0, :V23 => 0.0,
        :VL1 => 0.0, :VL2 => 0.0, :VL3 => 0.0
    )
end

function set_active!(settings::InputSettings, syms...)
    for s in syms
        push!(settings.enabled_inputs, s)
    end
end

function set_nominal!(settings::InputSettings, sym::Symbol, val::Float64)
    settings.nominal_values[sym] = val
end
