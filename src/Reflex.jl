"""
This file contains the definition of the reflex model, modified from both Heldt (2004) and Zamanian (2007). The model operates in continuous time to be compatible with the MTK architecture.
"""

"""
Afferent Arm
The afferent arm of the reflex model is defined by a linear state-space system with 9 states, acting as a tuned low-pass filter to account for the time delay in neuronal action. The output of the filter is compared with a set-point to find the error signal. An arctangent function is used to transform the error into a signal that can be saturated. The gain of the system is adjustable to account for differences between different afferent pathways.
"""

@component function Afferent(; name, A=A_9x9, B=B_9x1, C=C_1x9, p_set=p_cpr, gain=gain_cpr)
  @variables begin
      e(t)        # input
      y(t)        # output
      x(t)[1:9]   # state vector
      δ(t)        # Signal
  end

  @parameters begin
      p_set = p_set # Set Point (mmHg)
      gain = gain # Afferent Gain
  end

  D = Differential(t)
  eqs = []

  for i in 1:9
      push!(eqs, D(x[i]) ~ sum(A[i,j]*x[j] for j in 1:9) + B[i]*(e)) # Input Equation ẋ = Ax + Be
  end

  push!(eqs, y ~ sum(C[1,j]*x[j] for j in 1:9)) # Output Equation y = Cx (+ Du)
  push!(eqs, δ ~ gain*atan((y-p_set)/gain)) # Error signal saturation model

  ODESystem(eqs, t; name=name)
end

"""
Transfer Function
The following components build up the transfer function.
"""

"""
Pade Delay Approximation
The PadeDelay component implements a Pade approximation of a delay of length τ. Discrete delays are not easily implementable in the MTK architecture without introducing Delay Differential Equations, which cannot be structurally simplified. The order of the approximation is defined by the parameter n, which can be set externally. The delay is implemented as a state-space system with n states, where the output is the last state. Increasing n dynamically adjusts the initial conditions vector, but dramatically increases the computational time.
"""

@mtkmodel PadeDelay begin
  @structural_parameters begin
    n = 5
  end
  @parameters begin
    τ = 1.0
  end
  @variables begin
    u(t)
    y(t)
    x(t)[1:n]
  end
  @equations begin
    D(x[1]) ~ (u - x[1]) * (n/τ)
    [D(x[i]) ~ (x[i-1] - x[i]) * (n/τ) for i in 2:n]...
    y ~ x[n]
  end
end

"""
Double Integrator with Gain
The DoubleIntegratorGain component implements a double integrator with a gain. This is equivalent to a 1/s^2 transfer function, whose impulse response is a ramp function. The gain controls the slope of the ramp.
"""

@mtkmodel DoubleIntegratorGain begin
  @parameters begin
    gain = 1.0
  end
  @variables begin
    u(t)
    y(t)
    v(t)  # intermediate output
  end
  @equations begin
    D(v) ~ gain * u
    D(y) ~ v
  end
end

"""
Transfer Function Element
The TFElement component joins together the PadeDelay and DoubleIntegratorGain components to create a delay followed by a ramp.
"""

@mtkmodel TFElement begin
  @parameters begin
    τ = 1.0
    gain = 1.0
  end
  @structural_parameters begin
    delay_order = 5
  end
  @variables begin
    u(t)
    y(t)
  end
  @components begin
    tftime = PadeDelay(τ=τ, n=delay_order)
    double_integrator = DoubleIntegratorGain(gain=gain)
  end
  @equations begin
    u ~ tftime.u
    tftime.y ~ double_integrator.u
    y ~ double_integrator.y
  end
end

"""
Reflex Transfer Function
The TransferFunction component implements a transfer function that models the reflex response. Heldt (2004) models the reflex through its impulse response defined by a delay, peak, and end, scaled to have unit area. We mimic Zamanian by combining three delay + ramp elements. In sequence, the first delay + ramp begins to ramp up after tfdelay, the second delay + ramp then reverses this and begins the falling edge after time tfpeak. The final delay + ramp zeros out the system at time tfend. The gains are scaled dynamically to provide unit area, the gain equations are shown in the comments. A separate Impulse.jl file is included to test the impulse response of the system and the results of different Pade Approximation orders.
"""

@mtkmodel TransferFunction begin
  @structural_parameters begin
    delay_order = 5
  end
  @parameters begin
    reflex_delay = 2.0
    reflex_peak = 5.0
    reflex_end = 30.0
  end
  # g1 = 2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))
  # g2 = (-(2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) * (reflex_end - reflex_delay) / (reflex_end - reflex_peak))
  # g3 = -((2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) + (-(2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) * (reflex_end - reflex_delay) / (reflex_end - reflex_peak)))

  @components begin
    tfdelay = TFElement(τ=reflex_delay, gain=(2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))), delay_order=delay_order)
    tfpeak = TFElement(τ=reflex_peak, gain=-(2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) * (reflex_end - reflex_delay) / (reflex_end - reflex_peak), delay_order=delay_order)
    tfend = TFElement(τ=reflex_end, gain=-((2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) + (-(2 / ((reflex_end - reflex_delay) * (reflex_peak - reflex_delay))) * (reflex_end - reflex_delay) / (reflex_end - reflex_peak))), delay_order=delay_order)
  end
  @variables begin
    u(t)
    y(t)
  end
  @equations begin
    tfdelay.u ~ u
    tfpeak.u ~ u
    tfend.u ~ u
    y ~ tfdelay.y + tfpeak.y + tfend.y
  end
end

"""
NEW Reflexes
Here are the reflexes that are not included in the original model. These are the ones that are not included in the original model, but are necessary for the simulation.
"""

"""
Ursino: Autoregulation
"""

@mtkmodel CerebralAutoregulation begin
  @parameters begin
    _τO₂ = τO₂
    _gbO₂ = gbO₂
    _CvbO₂n = CvbO₂n
    _A_auto = A_auto
    _B_auto = B_auto
    _C_auto = C_auto
    _D_auto = D_auto
    _PaCO₂n = PaCO₂n
    _τCO₂ = τCO₂
  end
  @variables begin
    uCvbO₂(t)
    xbO₂(t)
    Phib(t)
    uPaCO₂(t)
    xbCO₂(t)
  end
  @equations begin
    D(xbO₂) ~ (-xbO₂ - _gbO₂ * (uCvbO₂ - _CvbO₂n)) / _τO₂
    Phib ~ ((_A_auto + (_B_auto) / (1 + _C_auto * exp(_D_auto * log(uPaCO₂)))) / (_A_auto + (_B_auto) / (1 + _C_auto * exp(_D_auto * log(_PaCO₂n))))) - 1
    D(xbCO₂) ~ (-xbCO₂ + Phib) / _τCO₂
  end
end

@mtkmodel Autoregulation begin
  @parameters begin
    _gjO₂ = 1.0
    _CvjO₂n = 1.0
    _τO₂ = 1.0
    _PaCO₂n = 1.0
    _kjCO₂ = 1.0
    _τCO₂ = 1.0
  end
  @variables begin
    uCvjO₂(t)
    xjO₂(t)
    uPaCO₂(t)
    Phij(t)
    xjCO₂(t)
  end
  @equations begin
    D(xjO₂) ~ (-xjO₂ - _gjO₂ * (uCvjO₂ - _CvjO₂n)) / _τO₂
    Phij ~ (1 - exp((uPaCO₂ - _PaCO₂n)/(_kjCO₂))) / (1 + exp((uPaCO₂ - _PaCO₂n)/(_kjCO₂)))
    D(xjCO₂) ~ (-xjCO₂ + Phij) / _τCO₂
  end
end

"""
Ursino: Chemoreceptors
"""

@mtkmodel Chemoreceptors begin
  @parameters begin
    Delay = 1.0
    Gain_A = 1.0
    Gain_f = 1.0
    set_point = 1.0
    time_A = 1.0
    time_f = 1.0
  end
  @structural_parameters begin
    delay_order = 5
  end
  @variables begin
    u(t)
    y_A(t)
    y_f(t)
  end
  @components begin
    delay = PadeDelay(τ=Delay, n=delay_order)
  end
  @equations begin
    u ~ delay.u
    D(y_A) ~ (-y_A + Gain_A * (delay.y - set_point)) / time_A
    D(y_f) ~ (-y_f + Gain_f * (delay.y - set_point)) / time_f
  end
end

"""
Ursino: Peripheral Chemoreceptors
"""

@mtkmodel PeripheralChemoreceptors begin
  @parameters begin
    _Ap = Ap
    _Bp = Bp
    _KO₂ = KO₂
    _KCO₂ = KCO₂ # (/s)
    _Cₜ = Cₜ # (ml/ml)
    _Kstat = Kstat # (/s)
    _τ_ph = τ_ph # (s)
    _τ_zh = τ_zh # (s)
    _τ_pl = τ_pl # (s)
    _Kdyn = Kdyn # (/s)
  end
  @variables begin
    uSaO₂(t)
    ucaCO₂(t)
    xO₂(t)
    ϕO₂(t)
    ϕCO₂(t)
    Phi(t)
    ϕCO₂dyn(t)
    ϕbarc(t)
    ϕc(t)
    fapc(t)
  end
  @equations begin
    xO₂ ~ _Ap * (1 - uSaO₂) + _Bp
    ϕO₂ ~ _KO₂ * (1 - exp(-xO₂ / _KO₂))
    ϕCO₂ ~ _KCO₂ * (ucaCO₂ - _Cₜ)
    Phi ~ ϕO₂ * ϕCO₂
    D(ϕCO₂dyn) ~ (_τ_zh * D(ucaCO₂) - ϕCO₂dyn) / _τ_ph
    ϕbarc ~ _Kstat * (1 - exp(-Phi / _Kstat)) + _Kdyn * (1 - exp(-ϕCO₂dyn / _Kdyn))
    D(ϕc) ~ (ϕbarc - ϕc) / _τ_pl
    fapc ~ ifelse(ϕc > 0, ϕc, 0.0)
  end
end

"""
Lung Stretch Receptors
"""

@mtkmodel LungStretch begin
  @parameters begin
    _Gasr = Gasr
    _τasr = τasr
  end
  @variables begin
    ϕasr(t)
    VT(t)
    fasr(t)
  end
  @equations begin
    ϕasr ~ _Gasr * VT
    D(fasr) ~ (-fasr + ϕasr) / _τasr
  end
end

"""
Ursino: Afferent Baroreflex
"""

@mtkmodel AfferentBaroreflex begin
  @parameters begin
    _τzb = τzb
    _τpb = τpb
    _Pn = Pn
    _kab = kab
    _fabₘₐₓ = fabₘₐₓ
    _fabₘᵢₙ = fabₘᵢₙ
  end
  @variables begin
    P(t)
    pb(t)
    fab(t)
  end
  @equations begin
    D(P) ~ (pb + (_τzb * D(pb)) - P) / _τpb
    fab ~ (_fabₘᵢₙ + _fabₘₐₓ * exp((P - _Pn) / _kab)) / (1 + exp((P - _Pn) / _kab))
  end
end

"""
Whittle: Afferent Cardiopulmonary Reflex
"""

@mtkmodel AfferentCPR begin
  @parameters begin
    _τzr = τzr
    _τpr = τpr
    _Prn = Prn
    _kcpr = kcpr
    _fcprₘₐₓ = fcprₘₐₓ
    _fcprₘᵢₙ = fcprₘᵢₙ
  end
  @variables begin
    P(t)
    pr(t)
    fcpr(t)
  end
  @equations begin
    D(P) ~ (pr + (_τzr * D(pr)) - P) / _τpr
    fcpr ~ (_fcprₘᵢₙ + _fcprₘₐₓ * exp((P - _Prn) / _kcpr)) / (1 + exp((P - _Prn) / _kcpr))
  end
end

"""
Ischemic Response
"""

@mtkmodel IschemicResponse begin
  @parameters begin
    _χₛⱼ = 1.0
    _PaO₂ₛⱼn = 1.0
    _kiscₛⱼ = 1.0
    _τisc = 1.0
    _PaCO₂n = 1.0
    _gccₛⱼ = 1.0
    _τcc = 1.0
    _θₛⱼn = 1.0
  end
  @variables begin
    uPaO₂(t)
    uPaCO₂(t)
    ωₛⱼ(t)
    ΔΘO₂ₛⱼ(t)
    ΔΘCO₂ₛⱼ(t)
    θₛⱼ(t)
  end
  @equations begin
    ωₛⱼ ~ _χₛⱼ / (1 + exp((uPaO₂ - _PaO₂ₛⱼn) / _kiscₛⱼ))
    D(ΔΘO₂ₛⱼ) ~ (-ΔΘO₂ₛⱼ + ωₛⱼ) / _τisc
    D(ΔΘCO₂ₛⱼ) ~ (-ΔΘCO₂ₛⱼ + _gccₛⱼ * (uPaCO₂ - _PaCO₂n)) / _τcc
    θₛⱼ ~ _θₛⱼn - ΔΘO₂ₛⱼ - ΔΘCO₂ₛⱼ
  end
end

"""
Ursino: Efferent Pathways
"""

@mtkmodel EfferentPathways begin
  @parameters begin
    # TODO
    _fes∞ = τ
    _fes₀ = gain
    _kes
    _Wb
    _Wc
    _Wp
  end
  @variables begin
    fab(t)
    fapc(t)
    fasr(t)
    fcpr(t)
    θₛₕ(t)
    θₛₚ(t)
    θₛᵥ(t)
    θᵥ(t)
    fₛₕ(t)
    fₛₚ(t)
    fₛᵥ(t)
    fᵥ(t)
  end
  @equations begin
    fₛₕ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbh * fab) + (_Wch * fapc) + (_Wph * fasr) - θₛₕ)), _fesₘₐₓ)
    fₛₚ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbp * fab) + (_Wcp * fapc) + (_Wpp * fasr) + (_Wrp * fcpr) - θₛₚ)), _fesₘₐₓ)
    fₛᵥ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbv * fab) + (_Wcv * fapc) + (_Wpv * fasr) + (_Wrv * fcpr) - θₛᵥ)), _fesₘₐₓ)

    fᵥ ~ (_fev₀ + _fev∞ * exp((fab - _fab₀) / _kev)) / (1 + exp((fab - _fab₀) / _kev)) + (_Wc_vagal * fapc) + (_Wp_vagal * fasr) - θᵥ
  end
end


