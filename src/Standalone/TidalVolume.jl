
using DifferentialEquations
using ModelingToolkit
using LinearAlgebra
using Symbolics
using Plots


@parameters t

D = Differential(t)

@connector PresPin begin
  p(t)
  q(t), [connect = Flow]
end

@mtkmodel OnePortPres begin
  @components begin
    out = PresPin()
    in = PresPin()
  end
  @variables begin
    Δp(t)
    q(t)
  end
  @equations begin
    Δp ~ out.p - in.p
    0 ~ in.q + out.q
    q ~ in.q
  end
end

cmH2O2mmHg = 0.73555912 # Conversion factor from cmH2O to mmHg

p₀ = 0.0 # Atmospheric pressure (mmHg)
pₚₗ₀ = -5.0 * cmH2O2mmHg # Initial Pleural pressure (mmHg)


"""
Lung Model
"""

#### Respiratory Muscle Parameters

p_musmin = -5.0 * cmH2O2mmHg # Minimum respiratory muscle pressure (mmHg)
RespRateₙₒₘ = 12.0 # Breathing Rate (breaths/min)
IE_ratio = 0.6 # Inspiratory to Expiratory Ratio

#### Lung Parameters

Rml = 1.021 * cmH2O2mmHg /1000 # Mouth-Larynx resistance (mmHg.s/ml: PRU)
Cl = 1.27 / cmH2O2mmHg # Larynx compliance (ml/mmHg)
V₀l = 34.4 # Larynx zero pressure volume (ml)
Rlt = 0.3369 * cmH2O2mmHg /1000 # Larynx-Trachea resistance (mmHg.s/ml: PRU)
Ctr = 2.38 / cmH2O2mmHg # Trachea compliance (ml/mmHg)
V₀tr = 6.63 # Trachea zero pressure volume (ml)
Rtb = 0.3063 * cmH2O2mmHg /1000 # Trachea-Bronchial resistance (mmHg.s/ml: PRU)
Cb = 13.1 / cmH2O2mmHg # Bronchial compliance (ml/mmHg)
V₀b = 18.7 # Brochea zero pressure volume (ml)
RbA = 0.0817 * cmH2O2mmHg /1000 # Bronchea-Alveolar resistance (mmHg.s/ml: PRU)
CA = 200 / cmH2O2mmHg # Alveolar compliance (ml/mmHg)
Ccw = 244.5 / cmH2O2mmHg # Chest wall compliance (ml/mmHg)
FRC = 2400 # Functional Residual Capacity (ml)

V₀A = FRC + (CA * pₚₗ₀) - (V₀l) - (Ctr * (p₀ - pₚₗ₀) + V₀tr) - (Cb * (p₀ - pₚₗ₀) + V₀b) # Alveolar

@mtkmodel ExternalPressureUB begin
  @components begin
    pext = PresPin()
  end
  @parameters begin
    p_ext = 0.0 # Baseline External pressure (mmHg)
  end
  @equations begin
    pext.p ~ p_ext
  end
end

@mtkmodel RespiratoryMuscles begin
  @extend OnePortPres()
  @parameters begin
    p = p_musmin
    RespRate₀ = RespRateₙₒₘ
    IEratio = IE_ratio
    ε = 0.0002         # width of the update pulse
    τ_reset = 0.01    # smoothing timescale
  end
  @variables begin
    RespRate_new(t) # breaths/min
    RespRate_chemo(t) # Chemoreceptor contributions (breaths/min)
    BreathInt_new(t) # Breath interval (s)
    ϕ(t) # Continuous slope
    BreathInt_held(t) # Breath interval held to a new breath (s)
    ϕ_wrapped(t) # Wrapped phase of the breathing cycle
    breath_trigger(t) # Trigger for the breath
    breath_trigger2(t) # Trigger for the end_inspiration

    t0(t) # Time withing the breathing cycle
    TE(t) # Expiration time (s)
    TI(t) # Inspiration time (s)
    τ(t) # Muscle time constant (s)
    p_new(t)
    p_chemo(t)
    p_held(t)
  end
  @equations begin
    RespRate_new ~ RespRate₀ + RespRate_chemo
    BreathInt_new ~ 60 / RespRate_new
    D(ϕ) ~ 1 / BreathInt_held
    ϕ_wrapped ~ ϕ - floor(ϕ)

    breath_trigger ~ exp(-((ϕ_wrapped - 1)^2) / ε)
    D(BreathInt_held) ~ breath_trigger * (BreathInt_new - BreathInt_held) / τ_reset

    t0 ~ ϕ_wrapped * BreathInt_held
    TE ~ BreathInt_held/(1 + IEratio)
    TI ~ TE * IEratio
    τ ~ TE/5

    p_new ~ p + p_chemo
    D(p_held) ~ breath_trigger * (p_new - p_held) / τ_reset

    breath_trigger2 ~ exp(-((ϕ_wrapped - TI / BreathInt_held)^2) / ε)


    # Time within the breathing cycle
    Δp ~ ifelse(t0 <= TI, (-p_held/(TI*TE)*t0^2 + p_held*BreathInt_held/(TI*TE)*t0),
          (p_held/(1-exp(-TE/τ))*(exp(-(t0-TI)/τ)-exp(-TE/τ))))
  end
end

"""
Lung
This is a complete model of the lung mechanics, based on the work of Albanese (2016). It has been modified to include the effects of altered-gravity and tilt on intrathoracic pressure.
"""

@mtkmodel Lung begin
  @components begin
    in = PresPin()
    chestwall = PresPin()
  end
  @parameters begin
    R_ml = Rml # Mouth-Larynx resistance (mmHg.s/ml: PRU)
    C_l = Cl # Larynx compliance (ml/mmHg)
    V₀_l = V₀l # Larynx zero pressure volume (ml)
    R_lt = Rlt # Larynx-Trachea resistance (mmHg.s/ml: PRU)
    C_tr = Ctr # Trachea compliance (ml/mmHg)
    V₀_tr = V₀tr # Trachea zero pressure volume (ml)
    R_tb = Rtb # Trachea-Bronchial resistance (mmHg.s/ml: PRU)
    C_b = Cb # Bronchial compliance (ml/mmHg)
    V₀_b = V₀b # Brochea zero pressure volume (ml)
    R_bA = RbA # Bronchea-Alveolar resistance (mmHg.s/ml: PRU)
    C_A = CA # Alveolar compliance (ml/mmHg)
    C_cw = Ccw # Chest wall compliance (ml/mmHg)
    V₀_A = V₀A # Alveolar zero pressure volume (ml)
    threshold = 1e-5
  end
  @variables begin
    pₐₒ(t) # External pressure (mmHg)
    p_l(t) # Larynx pressure (mmHg)
    p_tr(t) # Trachea pressure (mmHg)
    p_b(t) # Bronchial pressure (mmHg)
    p_A(t) # Alveolar pressure (mmHg)
    pₚₗ(t) # Pleural pressure (mmHg)
    V_l(t) # Larynx volume (ml)
    V_tr(t) # Trachea volume (ml)
    V_b(t) # Brochea volume (ml)
    V_A(t) # Alveolar volume (ml)
    V_D(t) # Dead space volume (ml)
    pₘᵤₛ(t) # Muscular pressure (mmHg)
    Vrᵢₙ(t) # Flow rate (ml/s)
    Vr_A(t) # Alveolar flow rate (ml/s)
    α(t) # Angle (radians)
    g(t) # Gravity (m/s^2)
    V_L(t)
  end
  @equations begin
    in.p ~ pₐₒ
    in.q ~ Vrᵢₙ

    chestwall.q ~ 0
    chestwall.p ~ pₘᵤₛ

    C_l * D(p_l) ~ (pₐₒ - p_l) / R_ml - (p_l - p_tr) / R_lt
    C_tr * (D(p_tr) - D(pₚₗ)) ~ (p_l - p_tr) / R_lt - (p_tr - p_b) / R_tb
    C_b * (D(p_b) - D(pₚₗ)) ~ (p_tr - p_b) / R_tb - (p_b - p_A) / R_bA
    C_A * (D(p_A) - D(pₚₗ)) ~ (p_b - p_A) / R_bA
    C_cw * (D(pₚₗ) - D(pₘᵤₛ)) ~ (p_l - p_tr) / R_lt

    Vrᵢₙ ~ (pₐₒ - p_l) / R_ml
    Vr_A ~ (p_b - p_A) / R_bA

    V_l ~ C_l * p_l + V₀_l
    V_tr ~ C_tr * (p_tr - pₚₗ) + V₀_tr
    V_b ~ C_b * (p_b - pₚₗ) + V₀_b
    V_A ~ C_A * (p_A - pₚₗ) + V₀_A
    V_D ~ V_l + V_tr + V_b

    V_L ~ V_A + V_D
  end
end

@mtkmodel TidalVolume begin
  @parameters begin
    τ_reset = 0.01 # Smoothing timescale (s)
  end
  @variables begin
    V_L(t) # Lung volume (ml)
    VL_min(t) # Min lung volume (ml)
    VL_max(t) # Peak lung volume (ml)
    breath_trigger(t) # New pressure (mmHg)
    breath_trigger2(t) # New pressure (mmHg)
    VT(t) # Tidal volume (ml)
  end
  @equations begin
    D(VL_min) ~ breath_trigger * (V_L - VL_min) / τ_reset
    D(VL_max) ~ breath_trigger2 * (V_L - VL_max) / τ_reset
    VT ~ VL_max - VL_min
  end
end

@named External = ExternalPressureUB(p_ext=p₀)
@named RespMuscles = RespiratoryMuscles()
@named Lungs = Lung()
@named TV = TidalVolume()

circ_eqs = [
  connect(Lungs.chestwall, RespMuscles.out),
  connect(External.pext, Lungs.in, RespMuscles.in),
  RespMuscles.p_chemo ~ 0.1*t,
  RespMuscles.RespRate_chemo ~ 0.0,
  TV.breath_trigger ~ RespMuscles.breath_trigger,
  TV.breath_trigger2 ~ RespMuscles.breath_trigger2,
  TV.V_L ~ Lungs.V_L,
]

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [
  External, Lungs, RespMuscles, TV, # Lung Model Breathing ChestWall
  ])

circ_sys = structural_simplify(circ_model)

equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_model))

u0 = [
  Lungs.p_l => p₀,
  Lungs.p_tr => p₀,
  Lungs.p_b => p₀,
  Lungs.p_A => p₀,
  RespMuscles.ϕ => 0.0,
  Lungs.pₚₗ => pₚₗ₀,
  RespMuscles.BreathInt_held => 60 / RespRateₙₒₘ,
  RespMuscles.p_held => p_musmin,

  TV.VL_min => 2000,
  TV.VL_max => 2500,

]

Start_time = 0.0
Time_step = 0.01
Stop_time = 30.0
tspan = (Start_time, Stop_time)

prob = ODEProblem(circ_sys, u0, tspan)

@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

display(plot(Sol, idxs=[TV.VT]))
display(plot(Sol, idxs=[RespMuscles.p_new]))