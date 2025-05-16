"""
This file contains the definition of the reflex model, modified from both Heldt (2004) and Zamanian (2007). The model operates in continuous time to be compatible with the MTK architecture.
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
    _fes₀ = fes₀
    _fesₘₐₓ = fesₘₐₓ
    _fes∞ = fes∞
    _kes = kes

    _Wbh = Wbₛₕ
    _Wch = Wcₛₕ
    _Wph = Wpₛₕ
    _Wrh = Wrₛₕ

    _Wbp = Wbₛₚ
    _Wcp = Wcₛₚ
    _Wpp = Wpₛₚ
    _Wrp = Wrₛₚ

    _Wbv = Wbₛᵥ
    _Wcv = Wcₛᵥ
    _Wpv = Wpₛᵥ
    _Wrv = Wrₛᵥ

    _fev₀ = fev₀
    _fev∞ = fev∞
    _kev = kev

    _fab₀ = fab₀
    _Wc_vagal = Wcᵥ
    _Wp_vagal = Wpᵥ
    _θᵥ = θᵥ
  end
  @variables begin
    fab(t)
    fapc(t)
    fasr(t)
    fcpr(t)
    θₛₕ(t)
    θₛₚ(t)
    θₛᵥ(t)
    fₛₕ(t)
    fₛₚ(t)
    fₛᵥ(t)
    fᵥ(t)
  end
  @equations begin
    fₛₕ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbh * fab) + (_Wch * fapc) + (_Wph * fasr) + (_Wrh * fcpr) - θₛₕ)), _fesₘₐₓ)
    fₛₚ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbp * fab) + (_Wcp * fapc) + (_Wpp * fasr) + (_Wrp * fcpr) - θₛₚ)), _fesₘₐₓ)
    fₛᵥ ~ min(_fes∞ + (_fes₀ - _fes∞) * exp(_kes * ((_Wbv * fab) + (_Wcv * fapc) + (_Wpv * fasr) + (_Wrv * fcpr) - θₛᵥ)), _fesₘₐₓ)

    fᵥ ~ (_fev₀ + _fev∞ * exp((fab - _fab₀) / _kev)) / (1 + exp((fab - _fab₀) / _kev)) + (_Wc_vagal * fapc) + (_Wp_vagal * fasr) - _θᵥ
  end
end

"""
Effectors
"""

@mtkmodel Effectors begin
  @parameters begin
    Gain = 1.0
    delay = 1.0
    time = 1.0
    min = 0.0
  end
  @structural_parameters begin
    delay_order = 2
  end
  @variables begin
    u(t)
    σθ(t)
    Δσ(t)
  end
  @components begin
    d = PadeDelay(τ=delay, n=delay_order)
  end
  @equations begin
    u ~ d.u
    σθ ~ ifelse(d.y >= min, Gain * log(d.y - min + 1), 0)
    D(Δσ) ~ (-Δσ + σθ) / time
  end
end

@mtkmodel EffectorsRR begin
  @parameters begin
    Gainₛ = 1.0
    Gainᵥ = 1.0
    delayₛ = 1.0
    delayᵥ = 1.0
    timeₛ = 1.0
    timeᵥ = 1.0
    min = 0.0
  end
  @structural_parameters begin
    delay_order = 2
  end
  @variables begin
    uₛ(t)
    uᵥ(t)
    σTₛ(t)
    σTᵥ(t)
    ΔTₛ(t)
    ΔTᵥ(t)
    ΔT(t)
  end
  @components begin
    dₛ = PadeDelay(τ=delayₛ, n=delay_order)
    dᵥ = PadeDelay(τ=delayᵥ, n=delay_order)
  end
  @equations begin
    uₛ ~ dₛ.u
    σTₛ ~ ifelse(dₛ.y >= min, Gainₛ * log(dₛ.y - min + 1), 0)
    D(ΔTₛ) ~ (-ΔTₛ + σTₛ) / timeₛ
    uᵥ ~ dᵥ.u
    σTᵥ ~ Gainᵥ * dᵥ.y
    D(ΔTᵥ) ~ (-ΔTᵥ + σTᵥ) / timeᵥ
    ΔT ~ ΔTₛ + ΔTᵥ
  end
end