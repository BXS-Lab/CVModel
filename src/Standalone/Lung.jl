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

v0pa = 160.0 # Pulmonary Artery zero pressure volume (ml)


HRₙₒₘ = 67.0  # Nominal Heart Rate (bpm)
TBV = 5300.0 # Total Blood Volume (ml)
ρ_b = 1060.0 # Blood Density (kg/m^3)
Pa2mmHg = 0.00750062 # Conversion factor from Pascal to mmHg
mmHg2dynecm2 = 10/Pa2mmHg # Conversion factor from mmHg to dyne/cm^2
cmH2O2mmHg = 0.73555912 # Conversion factor from cmH2O to mmHg
η_b = 0.004 # Blood Viscosity (Pa.s)

p₀ = 0.0 # Atmospheric pressure (mmHg)
pₚₗ₀ = -5.0 * cmH2O2mmHg # Initial Pleural pressure (mmHg)"

#### Respiratory Muscle Parameters

p_musmin = -5.0 * cmH2O2mmHg # Minimum respiratory muscle pressure (mmHg)
RespRateₙₒₘ = 12.0 # Breathing Rate (breaths/min)
IEratio = 0.6 # Inspiratory to Expiratory Ratio
Tbreath = 60.0 / RespRateₙₒₘ # Total Breathing Cycle Time (s)
T_E = Tbreath/(1 + IEratio) # Expiratory Time (s)
T_I = T_E * IEratio # Inspiratory Time (s)
τ_mus = T_E/5 # Respiratory muscle time constant (s)

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

V₀A = FRC + (CA * pₚₗ₀) - (V₀l) - (Ctr * (p₀ - pₚₗ₀) + V₀tr) - (Cb * (p₀ - pₚₗ₀) + V₀b) # Alveolar zero pressure volume (ml)

#### Tissue Volumes
Vₜ_brain = 1300 # Brain tissue volume (ml)
Vₜ_heart = 284 # Heart tissue volume (ml)
Vₜ_splanchnic = 2673 # Lung tissue volume (ml)
Vₜ_renal = 262 # Renal tissue volume (ml)
Vₜ_sm = 31200 # Skeletal muscle tissue volume (ml)
ub_mass_frac = 0.429 # Fraction of skeletal muscle in the upper body (Jansen 2000, https://journals.physiology.org/doi/full/10.1152/jappl.2000.89.1.81)
Vₜ_ub = Vₜ_sm * ub_mass_frac # Upper body (skeletal muscle) tissue volume (ml)
Vₜ_legs = Vₜ_sm * (1-ub_mass_frac) # Leg (skeletal muscle) tissue volume (ml)

#### Tissue O₂ Consumption Rate
MO₂_total = 250 # Total O₂ consumption rate (ml/min)
MO₂_spes_rat = 7.384 # Splanchnic to extrasplanchnic O₂ consumption rate ratio (Albanese, 2016)

MO₂_brain = 47.502/60 # Brain O₂ consumption rate (ml/s)
MO₂_heart = 24/60 # Heart O₂ consumption rate (ml/s)
MO₂_sm = 51.6/60 # Skeletal muscle O₂ consumption rate (ml/s)
MO₂_ub = MO₂_sm * ub_mass_frac # Upper body (skeletal muscle) O₂ consumption rate (ml/s)
MO₂_legs = MO₂_sm * (1-ub_mass_frac) # Leg (skeletal muscle) O₂ consumption rate (ml/s)
MO₂_rem = MO₂_total/60 - (MO₂_brain + MO₂_heart + MO₂_ub + MO₂_legs) # Remaining O₂ consumption rate (ml/min)
MO₂_renal = MO₂_rem / (1 + MO₂_spes_rat) # Renal O₂ consumption rate (ml/s)
MO₂_splanchnic = MO₂_renal * MO₂_spes_rat # Splanchnic O₂ consumption rate (ml/s)

#### Tissue CO₂ Production Rate
RQ = 0.84 # Respiratory Quotient (CO₂ produced / O₂ consumed)
MCO₂_total = MO₂_total * RQ # Total CO₂ production rate (ml/min)

MCO₂_brain = MO₂_brain * RQ # Brain CO₂ production rate (ml/s)
MCO₂_heart = MO₂_heart * RQ # Heart CO₂ production rate (ml/s)
MCO₂_splanchnic = MO₂_splanchnic * RQ # Splanchnic CO₂ production rate (ml/s)
MCO₂_renal = MO₂_renal * RQ # Renal CO₂ production rate (ml/s)
MCO₂_ub = MO₂_ub * RQ # Upper body (skeletal muscle) CO₂ production rate (ml/s)
MCO₂_legs = MO₂_legs * RQ # Leg (skeletal muscle) CO₂ production rate (ml/s)

"""
Lung Gas Exchange Model
"""
#### Environmental Conditions
FIO₂ = 21.0379 # Fraction of O₂ in the inspired air (%)
FICO₂ = 0.0421 # Fraction of CO₂ in the inspired air (%)
K = 1.2103 # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
pₐₜₘ = 760.0 # Atmospheric pressure (mmHg)
p_ws = 47.0 # Water vapor pressure (mmHg)

#### O₂ Dissociation Curve
CₛₐₜO₂ = 9.0 # Saturation of O₂ in the lungs (mmol/l)
h₁ = 0.3836
α₁ = 0.03198 # (/mmHg)
β₁ = 0.008275 # (/mmHg)
K₁ = 14.99 # (mmHg)

#### CO₂ Dissociation Curve
CₛₐₜCO₂ = 86.11 # Saturation of CO₂ in the lungs (mmol/l)
h₂ = 1.819
α₂ = 0.05591 # (/mmHg)
β₂ = 0.03255 # (/mmHg)
K₂ = 194.4 # (mmHg)

#### Physiological Status
Hgb = 15.0 # Hemoglobin concentration (g/dl)
sh = 0.017 # Pulmonary shunt fraction (fraction of blood that bypasses the lungs)

#### Constants
Hgb_O₂_binding = 1.34 # O₂ capacity (ml O₂/g Hgb), 1.34
sol_O₂ = 0.003/100 # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)

"""
This file contains the definitions of the compartment models used in the simulation. The underlying physiological equations are based on the work of Heldt (2004), with modifications by Zamanian (2007), Diaz Artiles (2015), and Whittle (2023). The Julia ModelingToolkit (MTK) is used to define the models, and the equations are expressed in a form suitable for numerical simulation through the DifferentialEquations.jl package.
"""

"""
Preamble
"""

@parameters t

D = Differential(t)

"""
Pin
This is a simple pin model with two variables: pressure (p, mmHg) and blood flow (q, ml/s). The flow is connected in accordance with Kirchhoff's laws.
"""

@connector Pin begin
  p(t)
  q(t), [connect = Flow]
  cO₂(t), [connect = Stream]
  cCO₂(t), [connect = Stream]
end

@connector PresPin begin
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
    cO₂(t)
    cCO₂(t)
  end
  @equations begin
    Δp ~ out.p - in.p
    0 ~ in.q + out.q
    q ~ in.q
    cO₂ ~ in.cO₂
    in.cO₂ ~ instream(in.cO₂)
    0 ~ cO₂ - out.cO₂
    cCO₂ ~ in.cCO₂
    in.cCO₂ ~ instream(in.cCO₂)
    0 ~ cCO₂ - out.cCO₂
  end
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

"""
Respiratory Muscles
This model represents the action of the respiratory muscles during breathing.
"""

@mtkmodel RespiratoryMuscles begin
  @extend OnePortPres()
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
  end
end



@mtkmodel LungGasExchange begin
  @parameters begin
    _FIO₂ = FIO₂ # Fractional concentration of O₂ in the inspired air (%)
    _FICO₂ = FICO₂ # Fractional concentration of CO₂ in the inspired air (%)
    _sh = sh # Shunt fraction
    _v0pa = v0pa # Unstressed volume of the pulmonary arterial blood (ml)
    _K = K # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
    _CₛₐₜO₂ = CₛₐₜO₂ # O₂ saturation constant (ml/ml)
    _CₛₐₜCO₂ = CₛₐₜCO₂ # CO₂ saturation constant (ml/ml)
    _h₁ = h₁ # O₂ saturation exponent
    _h₂ = h₂ # CO₂ saturation exponent
    _K₁ = K₁ # O₂ saturation constant (mmHg)
    _K₂ = K₂ # CO₂ saturation constant (mmHg)
    _α₁ = α₁ # O₂ saturation constant (/mmHg)
    _α₂ = α₂ # CO₂ saturation constant (/mmHg)
    _β₁ = β₁ # O₂ saturation constant (/mmHg)
    _β₂ = β₂ # CO₂ saturation constant (/mmHg)
    _pₐₜₘ = pₐₜₘ # Absolute Atmospheric pressure (mmHg)
    _p_ws = p_ws # Water vapor pressure (mmHg)
    _sol_O₂ = sol_O₂ # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)
    _Hgb_O₂_binding = Hgb_O₂_binding # Hemoglobin O₂ binding constant (ml O₂/ml blood/mmHg)
    _Hgb = Hgb # Hemoglobin concentration (g/dl)
  end
  @variables begin

    #### Connected from Lung Model
    Vrᵢₙ(t) # Flow rate (ml/s), must be connected externally
    Vr_A(t) # Alveolar flow rate (ml/s), must be connected externally
    V_D(t) # Dead space volume (ml), must be connected externally
    V_A(t) # Alveolar volume (ml), must be connected externally

    #### Connected from Pulmonary Artery
    qpa(t) # Pulmonary blood flow (ml/s), must be connected externally
    Vpa(t) # Pulmonary arterial blood volume (ml), must be connected externally
    cvO₂(t) # O₂ concentration in the deoxygenated blood (ml/ml), must be connected externally
    cvCO₂(t) # CO₂ concentration in the deoxygenated blood (ml/ml), must be connected externally


    FDO₂(t) # Fractional concentration of O₂ in the dead space (%)
    FDCO₂(t) # Fractional concentration of CO₂ in the dead space (%)
    FAO₂(t) # Fractional concentration of O₂ in the alveolar space (%)
    FACO₂(t) # Fractional concentration of CO₂ in the alveolar space (%)

    qpp(t) # Pulmonary blood flow (ml/s)
    qps(t) # Pulmonary shunted blood flow (ml/s)
    Vpp(t) # Pulmonary peripheral blood volume (ml)

    cppO₂(t) # O₂ concentration in the pulmonary blood (ml/ml)
    cppCO₂(t) # CO₂ concentration in the pulmonary blood (ml/ml)
    XppO₂(t) # O₂ partial pressure
    XppCO₂(t) # CO₂ partial pressure
    pppO₂(t) # O₂ partial pressure in pulmonary blood (mmHg)
    pppCO₂(t) # CO₂ partial pressure in pulmonary blood (mmHg)

    p_AO₂(t) # O₂ partial pressure in the alveolar space (mmHg)
    p_ACO₂(t) # CO₂ partial pressure in the alveolar space (mmHg)

    caO₂(t) # O₂ concentration in the arterial blood (ml/ml) (output)
    caCO₂(t) # CO₂ concentration in the arterial blood (ml/ml) (output)

    # XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
    # XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
    # paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
    # paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

    # SaO₂(t) # O₂ saturation in the arterial blood (%)


  end
  @equations begin

    #### Conservation of mass equations
    D(FDO₂) ~ (0.5 * (1 + tanh(100 * Vrᵢₙ)) * Vrᵢₙ * (_FIO₂ - FDO₂) +
           0.5 * (1 - tanh(100 * Vrᵢₙ)) * Vr_A * (FDO₂ - FAO₂)) / V_D
    D(FDCO₂) ~ (0.5 * (1 + tanh(100 * Vrᵢₙ)) * Vrᵢₙ * (_FICO₂ - FDCO₂) +
           0.5 * (1 - tanh(100 * Vrᵢₙ)) * Vr_A * (FDCO₂ - FACO₂)) / V_D
    D(FAO₂) ~ (0.5 * (1 + tanh(100 * Vrᵢₙ)) * Vr_A * (FDO₂ - FAO₂) - _K * (qpa * (1 - _sh) * (cppO₂ - cvO₂) + Vpp * D(cppO₂))) / V_A
    D(FACO₂) ~ (0.5 * (1 + tanh(100 * Vrᵢₙ)) * Vr_A * (FDCO₂ - FACO₂) - _K * (qpa * (1 - _sh) * (cppCO₂ - cvCO₂) + Vpp * D(cppCO₂))) / V_A
    Vpp ~ (Vpa - v0pa) * (1 - sh) + v0pa


    #### Dissociation equations
    cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
    XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
    cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
    XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

    #### Instantaneous equilibrium equations
    p_AO₂ ~ pppO₂
    p_ACO₂ ~ pppCO₂

    #### Gas fraction to partial pressure relationships
    p_AO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
    p_ACO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

    #### Mixing between capillary and shunted blood
    caO₂ ~ (qpp * cppO₂ + qps * cvO₂) / (qpp + qps)
    caCO₂ ~ (qpp * cppCO₂ + qps * cvCO₂) / (qpp + qps)

    qpp ~ qpa * (1 - _sh)
    qps ~ qpa * _sh

    #### O₂ saturation in arterial blood
    # caO₂ ~ CₛₐₜO₂ * (XaO₂)^(1/h₁)/(1 + (XaO₂)^(1/h₁))
    # XaO₂ ~ paO₂ * (1 + β₁ * paCO₂) / (K₁ * (1 + α₁ * paCO₂))
    # caCO₂ ~ CₛₐₜCO₂ * (XaCO₂)^(1/h₂)/(1 + (XaCO₂)^(1/h₂))
    # XaCO₂ ~ paCO₂ * (1 + β₂ * paO₂) / (K₂ * (1 + α₂ * paO₂))
    # SaO₂ ~ (caO₂ - paO₂ * sol_O₂) / (Hgb * Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
  end
end















"""
Model
"""

@named External = ExternalPressureUB()

@named RespMuscles = RespiratoryMuscles()
@named Lungs = Lung()
@named LungGE = LungGasExchange()

circ_eqs = [
  connect(External.pext, Lungs.in, RespMuscles.in),

  connect(Lungs.chestwall, RespMuscles.out),

  LungGE.Vrᵢₙ ~ Lungs.Vrᵢₙ,
  LungGE.Vr_A ~ Lungs.Vr_A,
  LungGE.V_D ~ Lungs.V_D,
  LungGE.V_A ~ Lungs.V_A,
  LungGE.qpa ~ 100,
  LungGE.Vpa ~ 150,
  LungGE.cvO₂ ~ 5,
  LungGE.cvCO₂ ~ 10
  # LungGE.caO₂ ~ Pulm_art.C.caO₂,
  # LungGE.caCO₂ ~ Pulm_art.C.caCO₂,
]

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [
  External, Lungs, RespMuscles, LungGE # Lung Model Breathing ChestWall
  ])

circ_sys = structural_simplify(circ_model)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
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

  LungGE.FDO₂ => FIO₂,
  LungGE.FDCO₂ => FICO₂,
  LungGE.FAO₂ => FIO₂,
  LungGE.FACO₂ => FICO₂
  # D(LungGE.pppO₂) => 0.0,
  # D(LungGE.pppCO₂) => 0.0
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

# @mtkmodel LungGasExchange begin
#   @parameters begin
#     _FIO₂ = FIO₂ # Fractional concentration of O₂ in the inspired air (%)
#     _FICO₂ = FICO₂ # Fractional concentration of CO₂ in the inspired air (%)
#     _sh = sh # Shunt fraction
#     _v0pa = v0pa # Unstressed volume of the pulmonary arterial blood (ml)
#     _K = K # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
#     _CₛₐₜO₂ = CₛₐₜO₂ # O₂ saturation constant (ml/ml)
#     _CₛₐₜCO₂ = CₛₐₜCO₂ # CO₂ saturation constant (ml/ml)
#     _h₁ = h₁ # O₂ saturation exponent
#     _h₂ = h₂ # CO₂ saturation exponent
#     _K₁ = K₁ # O₂ saturation constant (mmHg)
#     _K₂ = K₂ # CO₂ saturation constant (mmHg)
#     _α₁ = α₁ # O₂ saturation constant (/mmHg)
#     _α₂ = α₂ # CO₂ saturation constant (/mmHg)
#     _β₁ = β₁ # O₂ saturation constant (/mmHg)
#     _β₂ = β₂ # CO₂ saturation constant (/mmHg)
#     _pₐₜₘ = pₐₜₘ # Absolute Atmospheric pressure (mmHg)
#     _p_ws = p_ws # Water vapor pressure (mmHg)
#     _sol_O₂ = sol_O₂ # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)
#     _Hgb_O₂_binding = Hgb_O₂_binding # Hemoglobin O₂ binding constant (ml O₂/ml blood/mmHg)
#     _Hgb = Hgb # Hemoglobin concentration (g/dl)
#   end
#   @variables begin

#     #### Connected from Lung Model
#     Vrᵢₙ(t) # Flow rate (ml/s), must be connected externally
#     Vr_A(t) # Alveolar flow rate (ml/s), must be connected externally
#     V_D(t) # Dead space volume (ml), must be connected externally
#     V_A(t) # Alveolar volume (ml), must be connected externally

#     #### Connected from Pulmonary Artery
#     qpa(t) # Pulmonary blood flow (ml/s), must be connected externally
#     Vpa(t) # Pulmonary arterial blood volume (ml), must be connected externally
#     cvO₂(t) # O₂ concentration in the deoxygenated blood (ml/ml), must be connected externally
#     cvCO₂(t) # CO₂ concentration in the deoxygenated blood (ml/ml), must be connected externally


#     FDO₂(t) # Fractional concentration of O₂ in the dead space (%)
#     FDCO₂(t) # Fractional concentration of CO₂ in the dead space (%)
#     FAO₂(t) # Fractional concentration of O₂ in the alveolar space (%)
#     FACO₂(t) # Fractional concentration of CO₂ in the alveolar space (%)

#     qpp(t) # Pulmonary blood flow (ml/s)
#     qps(t) # Pulmonary shunted blood flow (ml/s)
#     Vpp(t) # Pulmonary peripheral blood volume (ml)

#     cppO₂(t) # O₂ concentration in the pulmonary blood (ml/ml)
#     cppCO₂(t) # CO₂ concentration in the pulmonary blood (ml/ml)
#     XppO₂(t) # O₂ partial pressure
#     XppCO₂(t) # CO₂ partial pressure
#     pppO₂(t) # O₂ partial pressure in pulmonary blood (mmHg)
#     pppCO₂(t) # CO₂ partial pressure in pulmonary blood (mmHg)

#     p_AO₂(t) # O₂ partial pressure in the alveolar space (mmHg)
#     p_ACO₂(t) # CO₂ partial pressure in the alveolar space (mmHg)

#     caO₂(t) # O₂ concentration in the arterial blood (ml/ml) (output)
#     caCO₂(t) # CO₂ concentration in the arterial blood (ml/ml) (output)

#     # XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
#     # XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
#     # paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
#     # paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

#     # SaO₂(t) # O₂ saturation in the arterial blood (%)


#   end
#   @equations begin

#     #### Conservation of mass equations
#     D(FDO₂) ~ ifelse(Vrᵢₙ > 0, Vrᵢₙ * (_FIO₂ - FDO₂), Vr_A * (FDO₂ - FAO₂)) / V_D
#     D(FDCO₂) ~ ifelse(Vrᵢₙ > 0, Vrᵢₙ * (_FICO₂ - FDCO₂), Vr_A * (FDCO₂ - FACO₂)) / V_D
#     D(FAO₂) ~ (ifelse(Vrᵢₙ > 0, Vr_A * (FDO₂ - FAO₂), 0) - _K * (qpa * (1 - _sh) * (cppO₂ - cvO₂) + Vpp * D(cppO₂))) / V_A
#     D(FACO₂) ~ (ifelse(Vrᵢₙ > 0, Vr_A * (FDCO₂ - FACO₂), 0) - _K * (qpa * (1 - _sh) * (cppCO₂ - cvCO₂) + Vpp * D(cppCO₂))) / V_A
#     Vpp ~ (Vpa - v0pa) * (1 - sh) + v0pa


#     #### Dissociation equations
#     cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
#     XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
#     cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
#     XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

#     #### Instantaneous equilibrium equations
#     p_AO₂ ~ pppO₂
#     p_ACO₂ ~ pppCO₂

#     #### Gas fraction to partial pressure relationships
#     FAO₂ ~ p_AO₂ / (_pₐₜₘ - _p_ws)
#     FACO₂ ~ p_ACO₂ / (_pₐₜₘ - _p_ws)

#     #### Mixing between capillary and shunted blood
#     caO₂ ~ (qpp * cppO₂ + qps * cvO₂) / (qpp + qps)
#     caCO₂ ~ (qpp * cppCO₂ + qps * cvCO₂) / (qpp + qps)

#     qpp ~ qpa * (1 - _sh)
#     qps ~ qpa * _sh

#     #### O₂ saturation in arterial blood
#     # caO₂ ~ CₛₐₜO₂ * (XaO₂)^(1/h₁)/(1 + (XaO₂)^(1/h₁))
#     # XaO₂ ~ paO₂ * (1 + β₁ * paCO₂) / (K₁ * (1 + α₁ * paCO₂))
#     # caCO₂ ~ CₛₐₜCO₂ * (XaCO₂)^(1/h₂)/(1 + (XaCO₂)^(1/h₂))
#     # XaCO₂ ~ paCO₂ * (1 + β₂ * paO₂) / (K₂ * (1 + α₂ * paO₂))
#     # SaO₂ ~ (caO₂ - paO₂ * sol_O₂) / (Hgb * Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
#   end
# end