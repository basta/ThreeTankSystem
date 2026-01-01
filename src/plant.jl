using ModelingToolkit

function build_nonlinear_plant()
    @independent_variables t
    @variables h1(t) h2(t) h3(t)
    @variables Q13V1(t) Q23V2(t) Q13V13(t) Q23V23(t) QL1(t) QL2(t) QN3(t)
    @parameters V1 = 0 V2 = 0 V13 = 0 V23 = 0 VL1 = 0 VL2 = 0 VL3 = 0
    @parameters Q1 = 0 Q2 = 0

    @parameters A = 0.0154 az = 1 Sh = 2e-5 g = 9.81 hv = 0.3 hmax = 0.62 Qimax = 1e-4 Ts = 5
    @parameters SL1 = 2e-5 SL2 = 2e-5 SL3 = 2e-5 SL13 = 2e-5 SL23 = 2e-5 S1 = 2e-5 S2 = 2e-5
    D = Differential(t)

    soft_sign_sqrt = x -> x / (sqrt(abs(x) + 1e-3))
    smooth_max = (x, y) -> 0.5 * (x + y + sqrt((x - y)^2 + 1e-4))
    smooth_max = (x, y) -> max(x, y) #override

    eqs = [
        D(h1) ~ 1 / A * (Q1 - Q13V1 - Q13V13 - QL1)
        D(h2) ~ 1 / A * (Q2 - Q23V2 - Q23V23 - QL2)
        D(h3) ~ 1 / A * (Q13V1 + Q13V13 + Q23V2 + Q23V23 - QN3)
        QL1 ~ VL1 * az * SL1 * soft_sign_sqrt(2 * g * h1)
        QL2 ~ VL2 * az * SL2 * soft_sign_sqrt(2 * g * h2)
        QN3 ~ VL3 * az * SL3 * soft_sign_sqrt(2 * g * h3)
        Q13V13 ~ V13 * az * SL13 * soft_sign_sqrt(2 * g * (h1 - h3))
        Q23V23 ~ V23 * az * SL23 * soft_sign_sqrt(2 * g * (h2 - h3))
        Q13V1 ~ V1 * az * S1 * soft_sign_sqrt(2 * g * (smooth_max(hv, h1) - smooth_max(hv, h3)))
        Q23V2 ~ V2 * az * S2 * soft_sign_sqrt(2 * g * (smooth_max(hv, h2) - smooth_max(hv, h3)))
    ]

    @named sys = ODESystem(eqs, t)
    return sys
end

