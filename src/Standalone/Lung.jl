"""
STANDALONE FILE – NOT CONNECTED TO THE MAIN CVModel.jl
This file is a standalone simulation of the lung model based on Albanese (2016). This model will be integrated into the main CVModel.jl package in the future.
"""

"""
Preamble
"""

using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using Symbolics

Start_time = 0.0
Time_step = 0.01
Stop_time = 30.0
tspan = (Start_time, Stop_time)

"""
Model Parameters
"""
cmH2O2mmHg = 0.73555912 # Conversion factor from cmH2O to mmHg

p_musmin = -5.0 * cmH2O2mmHg # Minimum respiratory muscle pressure (mmHg)
RespRateₙₒₘ = 12.0 # Breathing Rate (breaths/min)
IEratio = 0.6 # Inspiratory to Expiratory Ratio
Tbreath = 60.0 / RespRateₙₒₘ # Total Breathing Cycle Time (s)
T_E = Tbreath/(1 + IEratio) # Expiratory Time (s)
T_I = T_E * IEratio # Inspiratory Time (s)
τ_mus = T_E/5 # Respiratory muscle time constant (s)

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

p₀ = 0.0 # Atmospheric pressure (mmHg)
pₜₕ₀ = -5.0 * cmH2O2mmHg # Initial Intrathoracic (pleural) pressure (mmHg)

V₀A = FRC + (CA * pₜₕ₀) - (V₀l) - (Ctr * (p₀ - pₜₕ₀) + V₀tr) - (Cb * (p₀ - pₜₕ₀) + V₀b) # Alveolar zero pressure volume (ml)

"""
Compartment Definitions
"""

@parameters t

D = Differential(t)

@connector Pin begin
  p(t)
  q(t), [connect = Flow]
end

@mtkmodel OnePort begin
  @components begin
    out = Pin()
    in = Pin()
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

@mtkmodel RespiratoryMuscles begin
  @extend OnePort()
  @parameters begin
    p = p_musmin
    TI = T_I
    TE = T_E
    T0 = Tbreath
    τ = τ_mus
  end
  @variables begin
    t0(t) # Time withing the breathing cycle
    ϕ(t) # Phase of the breathing cycle
  end
  @equations begin
    D(ϕ) ~ 1 / T0
    t0 ~ (ϕ - floor(ϕ)) * T0
    # Δp ~ ifelse(t0 <= TI, 0.01 * t0, 0.1 * t0)
    # Time within the breathing cycle
    Δp ~ ifelse(t0 <= TI, (-p/(TI*TE)*t0^2 + p*T0/(TI*TE)*t0),
          (p/(1-exp(-TE/τ))*(exp(-(t0-TI)/τ)-exp(-TE/τ))))
  end
end

@mtkmodel ExternalPressureUB begin
  @components begin
    pext = Pin()
  end
  @parameters begin
    p_ext = 0.0 # Baseline External pressure (mmHg)
  end
  @equations begin
    pext.p ~ p_ext
  end
end

@mtkmodel Lung begin
  @components begin
    in = Pin()
    chestwall = Pin()
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
  end
  @variables begin
    pₐₒ(t) # External pressure (mmHg)
    p_l(t) # Larynx pressure (mmHg)
    p_tr(t) # Trachea pressure (mmHg)
    p_b(t) # Bronchial pressure (mmHg)
    p_A(t) # Alveolar pressure (mmHg)
    pₜₕ(t) # Intrathoracic (pleural) pressure (mmHg)
    V_l(t) # Larynx volume (ml)
    V_tr(t) # Trachea volume (ml)
    V_b(t) # Brochea volume (ml)
    V_A(t) # Alveolar volume (ml)
    V_D(t) # Dead space volume (ml)
    pₘᵤₛ(t) # Muscular pressure (mmHg)
    Vrᵢₙ(t) # Flow rate (ml/s)
    Vr_A(t) # Alveolar flow rate (ml/s)
  end
  @equations begin
    in.p ~ pₐₒ
    in.q ~ Vrᵢₙ

    chestwall.q ~ 0
    chestwall.p ~ pₘᵤₛ

    C_l * D(p_l) ~ (pₐₒ - p_l) / R_ml - (p_l - p_tr) / R_lt
    C_tr * (D(p_tr) - D(pₜₕ)) ~ (p_l - p_tr) / R_lt - (p_tr - p_b) / R_tb
    C_b * (D(p_b) - D(pₜₕ)) ~ (p_tr - p_b) / R_tb - (p_b - p_A) / R_bA
    C_A * (D(p_A) - D(pₜₕ)) ~ (p_b - p_A) / R_bA
    C_cw * (D(pₜₕ) - D(pₘᵤₛ)) ~ (p_l - p_tr) / R_lt

    Vrᵢₙ ~ (pₐₒ - p_l) / R_ml
    Vr_A ~ (p_b - p_A) / R_bA

    V_l ~ C_l * p_l + V₀_l
    V_tr ~ C_tr * (p_tr - pₜₕ) + V₀_tr
    V_b ~ C_b * (p_b - pₜₕ) + V₀_b
    V_A ~ C_A * (p_A - pₜₕ) + V₀_A
    V_D ~ V_l + V_tr + V_b
  end
end

"""
Model
"""

@named External = ExternalPressureUB(p_ext=p₀)

@named RespMuscles = RespiratoryMuscles()
@named Lung = Lung()

circ_eqs = [
  connect(External.pext, Lung.in, RespMuscles.in),
  connect(Lung.chestwall, RespMuscles.out)
]

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [
  External, Lung, RespMuscles # Lung Model Breathing ChestWall
  ])

circ_sys = structural_simplify(circ_model)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_model))

u0 = [
  Lung.p_l => p₀,
  Lung.p_tr => p₀,
  Lung.p_b => p₀,
  Lung.p_A => p₀,
  RespMuscles.ϕ => 0.0,
  Lung.pₜₕ => pₜₕ₀
]

prob = ODEProblem(circ_sys, u0, tspan)

@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

p1 = plot(Sol, idxs=[RespMuscles.out.p], xlims = (0, 15),
        label = "pₘᵤₛ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Respiratory Muscle Pressure")

p2 = plot(Sol, idxs=[Lung.pₜₕ], xlims = (0, 15),
        label = "pₜₕ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Intrathoracic Pressure")

p3 = plot(Sol, idxs=[Lung.p_A], xlims = (0, 15),
        label = "p_A",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Alveolar Pressure")

p4 = plot(Sol, idxs=[Lung.Vrᵢₙ], xlims = (0, 15),
        label = "Air Flow",
        xlabel = "Time (s)",
        ylabel = "Flow (ml/s)",
        title = "Air Flow")

p5 = plot(Sol, idxs=[Lung.V_A + Lung.V_D], xlims = (0, 15),
        label = "V_L",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Lung Volume")

p6 = plot(Sol, idxs=[Lung.V_A], xlims = (0, 15),
        label = "V_A",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Alveolar Volume")

p7 = plot(Sol, idxs=[Lung.V_D], xlims = (0, 15),
        label = "V_D",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Dead Space Volume")

display(plot(p1,p2,p3,p4,p5,p6,p7, layout=(4,2), size=(900,600), suptitle="Lungs"))

