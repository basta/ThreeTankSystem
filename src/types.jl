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
