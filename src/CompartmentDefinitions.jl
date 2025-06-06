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
This is a simple pin model with two variables: pressure (p, mmHg) and blood flow (q, ml/s). The flow is connected in accordance with Kirchhoff's laws. O₂ and CO₂ concentrations are also included as stream variables, which are carried by the flow. A separate PresPin model is defined for pressure-only connections (e.g., in the lung mechanics or tissue pressures).
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

"""
Ground
Ground is a single pin model that serves as a reference point for pressure. It has a single parameter P (mmHg), which represents the pressure at the ground node.
"""

@mtkmodel Ground begin
  @components begin
    g = Pin()
  end
  @parameters begin
    P = 0.0
  end
  @equations begin
    g.p ~ P
  end
end

"""
One Port
The One Port is a basic circuit element with an input and output. It has two variables: pressure difference (Δp, mmHg) and flow (q, ml/s). The pressure difference is defined as the difference between the output and input pressures. The flow is conserved such that the sum of the input and output flows is zero.
"""

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

"""
One Port with External Pressure
This model extends the One Port model by adding an external pressure pin (ep). The (variable) external pressure is connected to the output pin, and no flow is allowed through the external pressure pin.
"""

@mtkmodel OnePortWithExtPressure begin
  @components begin
    out = Pin()
    in = Pin()
    ep = PresPin()
  end
  @variables begin
    Δp(t)
    q(t)
    pg(t)
    cO₂(t)
    cCO₂(t)
  end
  @equations begin
    Δp ~ out.p - in.p
    0 ~ in.q + out.q
    0 ~ ep.q
    q ~ in.q
    pg ~ p - ep.p
    cO₂ ~ in.cO₂
    in.cO₂ ~ instream(in.cO₂)
    0 ~ cO₂ - out.cO₂
    cCO₂ ~ in.cCO₂
    in.cCO₂ ~ instream(in.cCO₂)
    0 ~ cCO₂ - out.cCO₂
  end
end

"""
Resistor
The Resistor model extends the One Port model and represents a simple resistive element. It has a single static parameter R (mmHg*s/ml = PRU), which represents the resistance of the element. The flow is related to the pressure difference by Ohm's law, where the flow is proportional to the negative pressure difference across the resistor.
"""

@mtkmodel Resistor begin
  @extend OnePort()

  @parameters begin
    R = 1.0
  end
  @equations begin
      Δp ~ -q * R
  end
end

"""
Variable Resistor
The Variable Resistor model extends the One Port model and represents a resistive element with a time-varying resistance. The resistance (mmHg*s/ml = PRU) is defined as a function of time, R(t), and can be externally controlled. Note that the equation is reversed such that q is the dependent variable and Δp is the independent variable. This is done to avoid difficulties with the structural simplify during implementation.
"""

@mtkmodel VarResistor begin
  @extend OnePort()

  @variables begin
    R(t) # Time-varying resistance
  end

  @equations begin
    q ~ -Δp / R
  end
end

"""
Variable Conductor
Same as the variable resistor but in conductance (G=1/R). This is used in the head capillaries where cerebral autoregulation affects conductance in a linear fashion.
"""

@mtkmodel VarConductor begin
  @extend OnePort()

  @variables begin
    G(t) # Time-varying resistance
  end

  @equations begin
    q ~ -Δp * G
  end
end

"""
Vessel Collapse models
These are resistors that exhibit vessel collapse behavior as defined by Zamanian (2007). Separate models are required for veins and arteries with the pressure-flow relationship reversed in order to avoid algebraic loops.
"""

@mtkmodel VesselCollapseVein begin
  @extend OnePort()
  @parameters begin
    R₀ = 1e-6
    ε = 1e-5
    Rₕᵢ = 1000
  end
  @variables begin
    R(t) # Time-varying resistance
    V(t) # Time-varying volume
    V₀(t)
  end

  @equations begin
    R ~ R₀ + Rₕᵢ * exp(log(ε/Rₕᵢ)*(V/V₀))
    q ~ -Δp / R
  end
end

@mtkmodel VesselCollapseArtery begin
  @extend OnePort()
  @parameters begin
    R₀ = 1e-6
    ε = 1e-5
    Rₕᵢ = 1000
  end
  @variables begin
    R(t) # Time-varying resistance
    V(t) # Time-varying volume
    V₀(t)
  end

  @equations begin
    R ~ R₀ + Rₕᵢ * exp(log(ε/Rₕᵢ)*(V/V₀))
    Δp ~ -q * R
  end
end


"""
Resistor Diode
The Resistor Diode model extends the One Port model and represents a resistive element with a diode-like behavior (e.g., a valve). The resistance (mmHg*s/ml = PRU) is defined as a static parameter R, and the flow set to zero when the pressure difference is negative.
"""

@mtkmodel ResistorDiode begin
  @extend OnePort()
  @parameters begin
    R = 1e-3
  end
  @equations begin
    q ~ -Δp / R * (Δp < 0)
  end
end

"""
Vessel Collapse Diode
The Vessel Collapse Diode model extends the resistor diode model to include a vessel collapse mechanic based on volume in a connected compliance. It is used in the arteries and veins with valves.
"""

@mtkmodel VesselCollapseDiode begin
  @extend OnePort()
  @parameters begin
    R₀ = 1e-6
    ε = 1e-5
    Rₕᵢ = 1000
  end
  @variables begin
    R(t) # Time-varying resistance
    V(t) # Time-varying volume
    V₀(t)
  end

  @equations begin
    R ~ R₀ + Rₕᵢ * exp(log(ε/Rₕᵢ)*(V/V₀))
    q ~ -Δp / R * (Δp < 0)
  end
end

"""
Starling Lung Resistor
This model represents the blood flow through the lungs as a series of parallel Starling Resistors. Based on 20% of the lung parenchyma below the pulmonary artery and vein, the model defines 3 zones for blood flow: continuous flow, periodic flow, and no flow. Equations from Heldt (2004).
"""

@mtkmodel StarlingResistor begin
  @components begin
    out = Pin()
    in = Pin()
  end
  @parameters begin
    R = 1.0
    h = 1.0
    ρ = ρ_b
    # p₀ = pₐₗᵥ
    ϵ = 1e-6  # small positive value
    Pa2mmHg_conv = Pa2mmHg # conversion factor from pascals to mmHg
  end
  @variables begin
    q(t)
    safe_gsinα(t)
    l₁(t)
    l₂(t)
    α(t)
    g(t)
    pₐₗᵥ(t)
    cO₂(t)
    cCO₂(t)
  end
  @equations begin
    0 ~ in.q + out.q
    q ~ in.q
    safe_gsinα ~ ifelse(abs(g * sin(α)) < ϵ, ϵ, g * sin(α))
    l₁ ~ clamp((out.p - pₐₗᵥ) / (ρ * safe_gsinα * Pa2mmHg_conv), -(h/100)/5, 4*(h/100)/5)
    l₂ ~ clamp((in.p - pₐₗᵥ) / (ρ * safe_gsinα * Pa2mmHg_conv), -(h/100)/5, 4*(h/100)/5)
    q ~ ((in.p - out.p) / ((h/100) * R) * (l₁ + (h/100) / 5)) + ((in.p - pₐₗᵥ) / ((h/100) * R) * (l₂ - l₁)) - ((ρ * g * Pa2mmHg_conv) / (2 * (h/100) * R) * sin(α) * (l₂^2 - l₁^2))
    cO₂ ~ in.cO₂
    in.cO₂ ~ instream(in.cO₂)
    0 ~ cO₂ - out.cO₂
    cCO₂ ~ in.cCO₂
    in.cCO₂ ~ instream(in.cCO₂)
    0 ~ cCO₂ - out.cCO₂
  end
end

"""
Capacitor
The Capacitor model extends the One Port model and represents a capacitive element. It has a single static parameter C (ml/mmHg), which represents the compliance of the element. The flow is related to the pressure difference by a first-order differential equation, where the rate of change of pressure difference is proportional to the negative flow through the capacitor.
"""

@mtkmodel Capacitor begin
  @extend OnePort()
  @parameters begin
          C = 1.0
  end
  @equations begin
          D(Δp) ~ -q / C
  end
end

"""
Inductor
The Inductor model extends the One Port model and represents an inductive element. It has a single static parameter L (mmHg*s^2/ml), which represents the inertia of the fluid in the system. The flow is related to the pressure difference by a first-order differential equation, where the rate of change of flow is proportional to the negative pressure difference across the inductor.
"""

@mtkmodel Inductance begin
  @extend OnePort()
  @parameters begin
    L = 1.0 # Inertia of the fluid in mmHg*s^2/ml
  end
  @equations begin
    D(q) ~ -Δp / L
  end
end

"""
Compliance
The Compliance model represents a vascular compliance. It has two parameters: the zero pressure volume (V₀, ml) and the compliance (C, ml/mmHg). Multiple flags are available to extend the model with additional features.
The inP flag indicates whether the formulation of the equations should be in terms of pressure (p) or volume (V). Pressure is used in this simulation, but volume is included for completeness.
The has_ep flag indicates whether there is an external static pressure, and the has_variable_ep flag extends this with a third pin to a time varying external pressure (both can be used at once for a bias). The variable p_rel is defined based on these two flags and the static parameter p₀ (mmHg).
The is_nonlinear flag can be set to transform the pressure-volume relationship into the nonlinear form as defined by Heldt (2004). If this flag is set, the V_max parameter (ml) is used to define the maximum distending volume. The nonlinear flag also inserts the variable qint (ml/s) to represent the flow through from these compartments to and from the interstitial compartment.
The has_abr and has_cpr flags introduce extra variables Vabr and Vcpr, which represent the reflex control of tone. The effective zero-pressure volume, V₀eff (ml) is a variable set by the sum of these and the base V₀. If introduced, Vabr and Vcpr must be connected externally. The variable pₜₘ (mmHg) is the transmural pressure, defined as the difference between the internal pressure and the external pressure.

Note: due to complexity this is composed as a @component and not a @mtkmodel. It makes no difference to the user.
"""

@component function Compliance(; name, V₀=0.0, C=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0, is_nonlinear=false, Flow_div = 1/3, V_max=1.0, V_min=1e-8, has_reflex=false, has_gasexchange=false, Vₜ=0.0, MO₂=0.0, RQ=0.84, is_pulmonary=false, pcol=p_col, is_heart=false)
  @named in = Pin() # Input pin
  @named out = Pin() # Output pin

  if has_variable_ep
    @named ep = PresPin() # If has_variable_ep, adds an external pressure pin
  end

  sts = @variables begin
    V(t) # Time-varying volume (ml)
    p(t) # Time-varying pressure (mmHg)
    V₀eff(t) # Time-varying Effective zero-pressure volume (ml)
    Cvar(t)
    cO₂(t) # Time-varying oxygen concentration (ml/ml)
    cCO₂(t) # Time-varying carbon dioxide concentration (ml/ml)
  end

  ps = @parameters begin
    V₀ = V₀ # Static zero-pressure volume (ml)
    C = C # Compliance (ml/mmHg)
    V_max = V_max # Maximum distending volume for nonlinear compartments (ml)
    Flow_div = Flow_div # Flow division to interstitial compartment (nonlinear only)
    V_min = V_min # Minimum volume (ml)
  end

  D = Differential(t)

  eqs = [
    0 ~ in.p - out.p,
    p ~ in.p
  ]

  if has_variable_ep # This statement defines the external pressure based on the flags
    push!(sts, (@variables p_rel(t))[1])
    if has_ep
      push!(ps, (@parameters p₀ = p₀)[1])
    end
    append!(eqs, [
      p_rel ~ ep.p + p₀,
      ep.q ~ 0
    ])
  elseif has_ep
    push!(ps, (@parameters p₀ = p₀)[1])
    p_rel = p₀
  else
    p_rel = p₀
  end

  if has_reflex # This statement defines the effective zero-pressure volume based on the flags
    push!(sts, (@variables ΔV(t))[1])
    append!(eqs, [V₀eff ~ max(V₀ + ΔV,V_min)])
  else
    append!(eqs, [V₀eff ~ V₀])
  end

  push!(sts, (@variables Cneg(t))[1])
append!(eqs, [Cneg ~ V₀eff / pcol])

  if is_nonlinear # Here are the nonlinear versions of the equations (including Qint)
    push!(sts, (@variables qint(t))[1])
    if inP
      append!(eqs, [ # Nonlinear in pressure
        Cvar ~ C / (1 + ((π * C) / (2 * V_max))^2 * (p - p_rel)^2), # Effective instantaneous compliance (ml/mmHg)
        V ~ V₀eff + (2 * V_max / π) * atan((π * C / (2 * V_max)) * (p - p_rel)),
        D(p) ~ (in.q + out.q - qint * Flow_div - D(V₀eff)) * 1 / Cvar + D(p_rel)
      ])
    else
      append!(eqs, [ # Nonlinear in volume
        Cvar ~ C / (1 + ((π * C) / (2 * V_max))^2 * ((2 * V_max / (π * C)) * tan((π * (V - V₀eff)) / (2 * V_max)))^2),# Effective instantaneous compliance (ml/mmHg)
        p ~ (2 * V_max / (π * C)) * tan((π * (V - V₀eff)) / (2 * V_max)) + p_rel,
        D(V) ~ in.q + out.q - qint * Flow_div
      ])
    end
  else # Here are the linear versions of the equations
    if inP # Linear in pressure
      append!(eqs, [
        # Cvar ~ C,
        Cvar ~ ifelse((p - p_rel) > 0, C, Cneg),
        V ~ (p - p_rel) * Cvar + V₀eff,
        D(p) ~ (1 / Cvar) * (in.q + out.q - D(V₀eff)) + D(p_rel) - 1/Cvar^2 * D(Cvar) * (V - V₀eff)
      ])
    else # Linear in volume
      append!(eqs, [
        Cvar ~ ifelse((V - V₀eff) > 0, C, Cneg),
        p ~ (V - V₀eff) / Cvar + p_rel,
        D(V) ~ in.q + out.q
      ])
    end
  end

  push!(sts, (@variables pₜₘ(t))[1])
  append!(eqs, [pₜₘ ~ p - p_rel])

  if has_gasexchange
    push!(ps, (@parameters Vₜ = Vₜ)[1])
    push!(ps, (@parameters MO₂ = MO₂)[1])
    push!(ps, (@parameters RQ = RQ)[1])
    append!(eqs, [
      cO₂ ~ in.cO₂,
      in.cO₂ ~ instream(in.cO₂),
      D(out.cO₂) ~ (in.q * (cO₂ - out.cO₂) - MO₂) / (V + Vₜ),
      cCO₂ ~ in.cCO₂,
      in.cCO₂ ~ instream(in.cCO₂),
      D(out.cCO₂) ~ (in.q * (cCO₂ - out.cCO₂) + MO₂ * RQ) / (V + Vₜ),
    ])
  elseif is_heart
    push!(ps, (@parameters Vₜ = Vₜ)[1])
    push!(sts, (@variables MO₂dyn(t))[1])
    push!(ps, (@parameters RQ = RQ)[1])
    append!(eqs, [
      cO₂ ~ in.cO₂,
      in.cO₂ ~ instream(in.cO₂),
      D(out.cO₂) ~ (in.q * (cO₂ - out.cO₂) - MO₂dyn) / (V + Vₜ),
      cCO₂ ~ in.cCO₂,
      in.cCO₂ ~ instream(in.cCO₂),
      D(out.cCO₂) ~ (in.q * (cCO₂ - out.cCO₂) + MO₂dyn * RQ) / (V + Vₜ),
    ])
  elseif is_pulmonary
    push!(sts, (@variables caO₂(t))[1])
    push!(sts, (@variables caCO₂(t))[1])
    append!(eqs, [
      in.cO₂ ~ instream(in.cO₂),
      cO₂ ~ in.cO₂,
      in.cCO₂ ~ instream(in.cCO₂),
      cCO₂ ~ in.cCO₂,
      out.cO₂ ~ caO₂,
      out.cCO₂ ~ caCO₂,
    ])
  else
    append!(eqs, [
      cO₂ ~ in.cO₂,
      in.cO₂ ~ instream(in.cO₂),
      D(out.cO₂) ~ in.q * (cO₂ - out.cO₂) / V,
      cCO₂ ~ in.cCO₂,
      in.cCO₂ ~ instream(in.cCO₂),
      D(out.cCO₂) ~ in.q * (cCO₂ - out.cCO₂) / V,
    ])
  end


  if has_variable_ep
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, ep)
  else
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
  end
end

"""
Sino-Atrial Node
This model represents the sino-atrial node of the heart. In our implementation it is a pacemaker that generates the beat waveform, φ_wrapped(t) based on the baseline RR interval (RR₀) and the baroreflex input representing both β-sympathetic and parasympathetic control. The RR interval is only updated at the start of a new beat, and a check on refractory period is also included. The atrial signal is provided separately as it is offset prior to ventricular contraction. We also pass the ABR β-sympathetic control of ventricular contractility through the SA node in order to hold until a new beat. These signals mediate the end-systolic elastance in the ventricles
"""

@mtkmodel SANode begin
  @parameters begin
    RR₀ = 1.0         # baseline RR interval (s)
    τ_reset = 0.01    # smoothing timescale
    ε = 0.002         # width of the update pulse
    τₐᵥ = 0.12        # AV delay (s)
    τ_refract = 0.2
    k = 100.0
  end

  @variables begin
    RRabr(t)          # baroreflex input
    RR_new(t)         # updated RR interval (instantaneous)
    RR_held(t)        # effective RR interval used by system
    ϕ(t)              # cumulative phase
    ϕ_wrapped(t)      # wrapped phase mod 1
    ϕ_wrapped_atria(t) # wrapped phase mod 1 with atrial offset
    atrial_offset(t)  # atrial offset
    beat_trigger(t)   # smoothed pulse function
    refractory_ok(t)
    Eabr_rv(t)
    Eabr_rv_held(t)
    Eabr_lv(t)
    Eabr_lv_held(t)
  end

  @equations begin
    RR_new ~ RR₀ + RRabr
    D(ϕ) ~ 1 / RR_held
    ϕ_wrapped ~ ϕ - floor(ϕ)

    ϕ_wrapped_atria ~ ϕ + τₐᵥ / sqrt(RR_held) - floor(ϕ + τₐᵥ / sqrt(RR_held))

    refractory_ok ~ 1 / (1 + exp(-k * RR_held * (ϕ_wrapped - τ_refract)))
    beat_trigger ~ exp(-((ϕ_wrapped - 1)^2) / ε) * refractory_ok

    D(RR_held) ~ beat_trigger * (RR_new - RR_held) / τ_reset
    D(Eabr_rv_held) ~ beat_trigger * (Eabr_rv - Eabr_rv_held) / τ_reset
    D(Eabr_lv_held) ~ beat_trigger * (Eabr_lv - Eabr_lv_held) / τ_reset
  end
end

"""
Heldt Chamber
This model represents a single chamber of the heart. It has two parameters: the zero pressure volume (V₀, ml) and the elastance (E, mmHg/ml). The elastance is time-varying and is defined by a shape function that describes the contraction and relaxation of the chamber. The model also includes a transmural pressure variable (pₜₘ, mmHg) that can be directly interrogated, an external static pressure (p₀, mmHg), and a time-varying external pressure (ep). The equations are defined in terms of either pressure or volume, depending on the inP flag. The end-systolic elastance can be adjusted based on the ABR signal (Eabr_held) from the SA node if the has_abr flag is set. The end-systolic elastance is defined as Eₘₐₓeff, which is limited by Elimit to prevent excessive elastance with severe orthostatic stress.

Note: due to complexity this is composed as a @component and not a @mtkmodel. It makes no difference to the user.
"""

@component function HeldtChamber(; name, inP=false, V₀, p₀=0.0, Eₘᵢₙ, Eₘₐₓ, Elimit=100.0, τₑₛ, has_abr=false)

  @named in = Pin()
  @named out = Pin()
  @named ep = PresPin()

  sts = @variables begin
    V(t) # Volume (ml)
    p(t) # Pressure (mmHg)
    pₜₘ(t) # Transmural pressure (mmHg)
    ϕ(t) # Fraction of beat cycle
    τ(t) # Instantaneous RR Interval (s)
    tᵢ(t) # Time within beat (s)
    τes(t) # End Systole Time (s)
    shape(t) # shape function for the time varying elastance
    dshape(t) # derivative of the shape function
    E(t) # Time varying elastance (mmHg/ml)
    DE(t) # Derivative of E(t)
    p_rel(t) # Time varying external pressure (mmHg)
    Eabr_held(t) # ABR ventricular contractility
    Eₘₐₓeff(t) # Beat-to-Beat end-systolic elastance (mmHg/ml)
    cO₂(t)
    cCO₂(t)
  end

  ps = @parameters begin
    V₀ = V₀ # Zero pressure volume (ml)
    p₀ = p₀ # External bias pressure (mmHg)
    Eₘᵢₙ = Eₘᵢₙ # Diastolic elastance (mmHg/ml)
    Eₘₐₓ = Eₘₐₓ # Baseline end-systolic elastance (mmHg/ml)
    Elimit = Elimit # Elastance limit to prevent excessive elastance with severe orthostatic stress (mmHg/ml)
    τₑₛ = τₑₛ     # end systole (rising edge end) time constant
  end

  D = Differential(t)
  eqs = Equation[]

  append!(eqs, [
    p_rel ~ ep.p + p₀,
    ep.q ~ 0,
    0 ~ in.p - out.p,
    p ~ in.p,
    pₜₘ ~ p - p_rel, # Transmural pressure

    tᵢ ~ ϕ * τ, # Modulated time within contraction cycle
  ])


    if has_abr # If there is ABR control, set the adjusted end-systolic elastance based on the held ABR signal
      append!(eqs, [
        Eₘₐₓeff ~ min(Eₘₐₓ + Eabr_held, Elimit),
      ])
    else # Otherwise static end-systolic elastance
      append!(eqs, [
        Eₘₐₓeff ~ Eₘₐₓ,
      ])
    end

    #### Time-varying elastance
    append!(eqs, [
    τes ~ (τₑₛ * sqrt(τ)), # Time varying end-systolic time (s)

    shape ~ (tᵢ <= τes) * (1 - cos((tᵢ / τes) * pi)) / 2 +
    (tᵢ > τes) * (tᵢ <= (1.5 * τes)) * (1 + cos(((tᵢ - τes) / τes) * 2 * pi)) / 2 +
    (tᵢ > (1.5 * τes)) * 0, # Shape function for the time varying elastance

    dshape ~ (tᵢ <= τes) * pi * sin((tᵢ / τes) * pi) / 2 * ((D(tᵢ) * τes - tᵢ * D(τes)) / τes^2) +
    (tᵢ > τes) * (tᵢ <= (1.5 * τes)) * pi * -sin(((tᵢ - τes) / τes) * 2 * pi) * ((D(tᵢ) * τes - (tᵢ - τes) * D(τes)) / τes^2) +
    (tᵢ > (1.5 * τes)) * 0, # Derivative of the shape function

    E ~ Eₘᵢₙ + (Eₘₐₓeff - Eₘᵢₙ) * shape, # Time varying elastance (mmHg/ml)

    DE ~ D(Eₘₐₓeff) * shape + (Eₘₐₓeff - Eₘᵢₙ) * dshape # Derivative of E(t)
  ])

  if inP # If inP is true, the equations are defined in terms of pressure
    append!(eqs, [
      V ~ (p - p_rel) / E + V₀,
      D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE + D(p_rel),
    ])
  else # Otherwise, the equations are defined in terms of volume
    append!(eqs, [
      p ~ (V - V₀) * E + p_rel,
      D(V) ~ in.q + out.q,
    ])
  end
  append!(eqs, [
    cO₂ ~ in.cO₂,
    in.cO₂ ~ instream(in.cO₂),
    D(out.cO₂) ~ in.q * (cO₂ - out.cO₂) / V,
    cCO₂ ~ in.cCO₂,
    in.cCO₂ ~ instream(in.cCO₂),
    D(out.cCO₂) ~ in.q * (cCO₂ - out.cCO₂) / V,
  ])


  compose(ODESystem(eqs, t, sts, ps; name=name), in, out, ep)
end

"""
Hydrostatic Pressure
This model represents the hydrostatic pressure in a fluid column. It has a static parameter ρ (kg/m^3), which represents the density of the fluid, and a vertical length (cm). There are two time-varying inputs: the angle α (radians) and the gravity level g (m/s^2). These must be connected externally. The 'con' parameter makes the pressure drop half the vertical height of the column, but can be adjusted (e.g., 3.0 in the leg compartments). The numerical constant in the equation converts the output to mmHg.
"""

@mtkmodel HydrostaticPressure begin
  @extend OnePort()
  @variables begin
    α(t)  # time-varying angle input (radians, to be connected externally)
    g(t)  # time-varying gravity input (m/s^2, to be connected externally)
  end
  @parameters begin
    ρ = ρ_b # density of the blood (kg/m^3)
    h = 10.0 # vertical length (cm)
    con = con_default # conversion factor (e.g., 2.0 for half height)
    Pa2mmHg_conv = Pa2mmHg # conversion factor from pascals to mmHg
  end
  @equations begin
    Δp ~ ρ * g * (h / con / 100) * sin(α) * Pa2mmHg_conv
  end
end

"""
Tissue Pressure
This model represents the tissue pressure on a fluid column. It has a static parameter ρ (kg/m^3), which represents the density of fat free tissue, and a compartment radius, rad (cm). There are two time-varying inputs: the angle α (radians) and the gravity level g (m/s^2). These must be connected externally. The numerical constant in the equation converts the output to mmHg.
"""

@mtkmodel TissuePressure begin
  @extend OnePortPres()
  @variables begin
    α(t)  # time-varying angle input (radians, to be connected externally)
    g(t)  # time-varying gravity input (m/s^2, to be connected externally)
  end
  @parameters begin
    ρ = ρ_fft # density of fat free tissue (kg/m^3)
    rad = 10.0 # Radius of compartment (cm)
    Pa2mmHg_conv = Pa2mmHg # conversion factor from pascals to mmHg
  end
  @equations begin
    Δp ~ ρ * g * (rad/100) * cos(α) * Pa2mmHg_conv # mmHg conversion
  end
end

"""
Artery
This model represents an arterial compartment. It is a lumped compartment consisting of a resistor and a compliance. If optional hydrostatic pressure is included it is added in series with the resistor. If optional tissue pressure is included it is added as an external force on the compliance. These must be externally connected to gravity and angle. There is also an optional flag to replace the resistor with a resistor diode combo (used for the pulmonary and aortic valves). The transmural pressure can be directly interrogated, along with the pressure drops across the resistor, hydrostatic pressure, and tissue pressure.
"""

@mtkmodel Artery begin
  @structural_parameters begin
    R = 1.0
    C = 1.0
    p₀ = 0.0
    V₀ = 1.0
    ρ = ρ_b
    h = 10.0
    con = 2.0
    has_valve = false
    has_hydrostatic = true  # new flag to include/exclude Ph
    has_tissue = true
    has_inertia = true
    ρt = ρ_fft
    rad = 10.0
    L = 1.0
    has_gasexchange=false
    Vₜ=0.0
    MO₂=0.0
    RQ=RQ₀
    is_pulmonary=false
    is_heart=false
  end

  @variables begin
    Δp(t)
    q(t)
    α(t)         # input angle
    g(t)         # input gravity
    Δp_Ph(t)     # pressure drop across hydrostatic pressure (set to 0 if not used)
    Δp_Pt(t)      # pressure drop across tissue pressure (set to 0 if not used)
    Δp_R(t)      # pressure drop across resistor
    pₜₘ(t)     # transmural pressure
    cO₂(t)
    cCO₂(t)
  end

  @components begin
    in = Pin()
    out = Pin()
    ep = PresPin()

    if has_valve
      R = VesselCollapseDiode(R₀=R)
    else
      R = VesselCollapseArtery(R₀=R)
    end
    C = Compliance(V₀=V₀, C=C, inP=true, has_ep=true, has_variable_ep=true, p₀=p₀, is_nonlinear=false, has_gasexchange=has_gasexchange, Vₜ=Vₜ, MO₂=MO₂, RQ=RQ, is_pulmonary=is_pulmonary, is_heart=is_heart)
    if has_hydrostatic
      Ph = HydrostaticPressure(ρ=ρ, h=h, con=con)
    end
    if has_tissue
      Pt = TissuePressure(ρ=ρt, rad=rad)
    end
    if has_inertia
      L = Inductance(L=L)
    end
  end

  @equations begin
    Δp ~ out.p - in.p
    q ~ in.q
    in.cO₂ ~ cO₂
    in.cCO₂ ~ cCO₂
    if has_inertia
      connect(in, L.in)
      connect(L.out, R.in)
    else
      connect(in, R.in)
    end
    if has_hydrostatic
      connect(R.out, Ph.in)
      connect(Ph.out, C.in)
      Ph.α ~ α
      Ph.g ~ g
      Δp_Ph ~ Ph.Δp
    else
      connect(R.out, C.in)
      Δp_Ph ~ 0  # no pressure drop if Ph is disabled
    end
    connect(C.out, out)
    R.V ~ C.V
    R.V₀ ~ C.V₀eff
    if has_tissue
      connect(C.ep, Pt.out)
      connect(Pt.in, ep)
      Pt.α ~ α
      Pt.g ~ g
      Δp_Pt ~ Pt.Δp
    else
      connect(C.ep, ep)
      Δp_Pt ~ 0  # no pressure drop if tissue is disabled
    end
    pₜₘ ~ C.pₜₘ
    Δp_R ~ R.Δp
  end
end

"""
Vein
This model represents a venous compartment. It is a lumped compartment consisting of a compliance and a resistor. If optional hydrostatic pressure is included it is added in series with the compliance. If optional tissue pressure is included it is added as an external force on the compliance. These must be externally connected to gravity and angle. There is also an optional flag to replace the resistor with a resistor diode combo. The transmural pressure can be directly interrogated, along with the pressure drops across the resistor, hydrostatic pressure, and tissue pressure. The compliance is optionally nonlinear, and the effective zero-pressure volume can be adjusted with the ABR and CPR signals Vabr and Vcpr (connected externally). The flow division to the interstitial compartment is connected here in the nonlinear case and feeds into the compliance.
"""


@mtkmodel Vein begin
  @structural_parameters begin
    R = 1.0
    C = 1.0
    p₀ = 0.0
    V₀ = 1.0
    V_min = 1e-8
    ρ = ρ_b
    h = 10.0
    con = 2.0
    has_valve = false
    has_hydrostatic = true  # new flag to control inclusion of hydrostatic component
    is_nonlinear = false
    V_max = 1.0
    Flow_div = 1/3
    has_reflex = false
    has_tissue = true
    ρt = ρ_fft
    rad = 10.0
  end

  @variables begin
    Δp(t)
    q(t)
    α(t)
    g(t)
    Δp_Ph(t)
    Δp_Pt(t)      # pressure drop across tissue pressure (set to 0 if not used)
    Δp_R(t)
    pₜₘ(t)
    cO₂(t)
    cCO₂(t)
    if has_reflex
      ΔV(t)
    end
  end

  @components begin
    in = Pin()
    out = Pin()
    ep = PresPin()

    C = Compliance(V₀=V₀, C=C, inP=true, has_ep=true, has_variable_ep=true, p₀=p₀, is_nonlinear=is_nonlinear, V_max=V_max, V_min=V_min, Flow_div=Flow_div, has_reflex=has_reflex)

    if has_hydrostatic
      Ph = HydrostaticPressure(ρ=ρ, h=h, con=con)
    end

    if has_valve
      R = VesselCollapseDiode(R₀=R)
    else
      R = VesselCollapseVein(R₀=R)
    end
    if has_tissue
      Pt = TissuePressure(ρ=ρt, rad=rad)
    end

  end

  @equations begin
    Δp ~ out.p - in.p
    q ~ in.q
    cO₂ ~ in.cO₂
    cCO₂ ~ in.cCO₂
    connect(in, C.in)
    connect(C.out, has_hydrostatic ? Ph.in : R.in)

    if has_tissue
      connect(C.ep, Pt.out)
      connect(Pt.in, ep)
      Pt.α ~ α
      Pt.g ~ g
      Δp_Pt ~ Pt.Δp
    else
      connect(C.ep, ep)
      Δp_Pt ~ 0  # no pressure drop if tissue is disabled
    end

    if has_reflex
      C.ΔV ~ ΔV
    end

    if has_hydrostatic
      connect(Ph.out, R.in)
      Ph.α ~ α
      Ph.g ~ g
      Δp_Ph ~ Ph.Δp
    else
      Δp_Ph ~ 0  # no hydrostatic drop
    end

    connect(R.out, out)
    R.V ~ C.V
    R.V₀ ~ C.V₀eff
    Δp_R ~ R.Δp
    pₜₘ ~ C.pₜₘ
  end
end

"""
Vertebral Plexus
Currently not implemented.
"""

@mtkmodel VertebralPlexus begin
  @structural_parameters begin
    R = 1.0        # Zero pressure volume (ml)
    ρ = ρ_b
    h = 10.0
    con = 2.0
  end

  @variables begin
    α(t)            # time-varying angle input (radians, to be connected externally)
    g(t)            # time-varying gravity input (m/s^2, to be connected externally)
    Δp(t)           # pressure drop across the compartment
    Δp_Ph(t)     # pressure drop across hydrostatic pressure (set to 0 if not used)
    q(t)            # flow rate into the compartment
    cO₂(t)
    cCO₂(t)
  end

  @components begin
    in = Pin()
    out = Pin()
    R = ResistorDiode(R=R)
    Ph = HydrostaticPressure(ρ=ρ, h=h, con=con)
  end

  @equations begin
    Δp ~ out.p - in.p
    q ~ in.q
    connect(in, R.in)
    connect(R.out, Ph.in)
    connect(Ph.out, out)
    Ph.α ~ α
    Ph.g ~ g
    Δp_Ph ~ Ph.Δp
    cO₂ ~ in.cO₂
    in.cO₂ ~ instream(in.cO₂)
    0 ~ cO₂ - out.cO₂
    cCO₂ ~ in.cCO₂
    in.cCO₂ ~ instream(in.cCO₂)
    0 ~ cCO₂ - out.cCO₂
  end
end

"""
Interstitial Compartment
This model represents the interstitial compartment. It is a lumped compartment acting as a fluid reserve filled by an RC network. The time constant is defined by the τ parameter (s). The interstitial compartment is connected to the nonlinear venous compartments through a flow variable Qint (ml/s). The maximum instantaneous volume of the interstitial compartment is defined by the Vmax variable (ml), additive for tilt and LBNP and defined at a reference of 85° HUT and/or –70 mmHg. The RC network acts as an exponential approach to the maximal value.
"""

@mtkmodel InterstitialCompartment begin
  @parameters begin
    Vmtilt = 1.0        # Scaling factor for Vmax in tilt
    Vmlbnp = 1.0        # Scaling factor for Vmax in LBNP
    τ = 276.0       # Time constant
  end

  @variables begin
    α(t)
    g(t)
    p_lbnp(t)
    Vmax(t)
    Qint(t)
    Vint(t)
  end

  @equations begin
    # Compute Vmax from alpha
    Vmax ~ Vmtilt * sin(α) / sin(85 / 180 * π) * (g / 9.81) + Vmlbnp * -1 * p_lbnp / 70

    τ * D(Qint) + Qint ~ D(Vmax)

    D(Vint) ~ (Vmax - Vint) / τ

  end
end

"""
Intrathoracic Pressure
This model represents the intrathoracic pressure. It currently just connects the cardivascular model to the lung model. Hydrostatic effects are WIP.
"""

@mtkmodel IntrathoracicPressure begin
  @components begin
    pth = PresPin()
  end
  # @parameters begin
  #   pₜₕ = -4.0 # Baseline intrathoracic pressure (mmHg)
  # end
  @variables begin
    α(t)
    g(t)
    pₚₗ(t)
  end
  @equations begin
    pth.p ~ pₚₗ# - 3.5 * (g / 9.81) * sin(α)
  end
end

"""
IntraAbdominal Pressure
This model represents the intra-abdominal pressure. It is defined by a baseline pressure (p_abd) and currently has no time-varying effects.
"""

@mtkmodel IntraAbdominalPressure begin
  @components begin
    pabd = PresPin()
  end
  @parameters begin
    p_abd = 0.0 # Baseline IntraAbdominal pressure (mmHg)
  end
  @equations begin
    pabd.p ~ p_abd
  end
end

"""
External Pressure
This model represents the external pressure. It is defined by a baseline pressure (p_ext) and currently has no time-varying effects.
"""

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
Intracranial Pressure
This model represents the intracranial pressure. It is defined by a baseline pressure (p_icp) and currently has no time-varying effects.
"""

@mtkmodel IntracranialPressure begin
  @components begin
    picp = PresPin()
  end
  @parameters begin
    p_icp = 10.0 # Baseline External pressure (mmHg)
  end
  @equations begin
    picp.p ~ p_icp
  end
end

"""
External Pressure (Lower Body, including LBNP)
This model represents the external pressure on the legs. It is defined by a baseline external pressure (p_ext) and is connected externally to the LBNP driver.
"""

@mtkmodel ExternalPressureLB begin
  @components begin
    pext = PresPin()
  end
  @variables begin
    p_lbnp(t)
  end
  @parameters begin
    p_ext = 0.0 # Baseline External pressure (mmHg)
  end
  @equations begin
    pext.p ~ p_ext + p_lbnp
  end
end

"""
`MynardValve_Atrioventricular(; name, ρ, Mrg, Mst, Ann, Kvc, Kvo)`

Implements the Mynard description for flow across the atrioventricular valves, full description in [Mynard].
This valve description corresponds to the atrioventricular valves where interia is not considered.

Note: The minimum level of regurgitation has to be set to machine precision eps()

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)
p is scaled to ensure units are consistent throughout.

Named parameters:

name    name of the element
`ρ`     Blood density in kg/m^3
`Mrg`   Level of regurgitation exhibited by a valve in DN
`Mst`   Level of stenosis exhibited by a valve in DN
`Ann`   Annulus area in cm^2
`Kvc`   Valve closing rate coefficent in  cm^2/(dynes*s)
`Kvo`   Valve opening rate coefficent in cm^2/(dynes*s)

p is calculated in mmHg
q is calculated in cm^3/s (ml/s)
"""

@mtkmodel MynardValve_Atrioventricular begin
  @extend OnePort()
  @parameters begin
          ρ = ρ_b
          Mrg = 0.0
          Mst = 1.0
          Ann
          Kvc
          Kvo
  end
  @variables begin
          Aeff(t)
          ζ(t)
          B(t)
          Aeff_min(t)
          Aeff_max(t)
          L(t)
  end
  begin
    Δp = -mmHg2dynecm2 * Δp # Convert mmHg to dynes/cm^2 = 1 Barye = 1 g/(cm*s^2)
  end
  @equations begin
          # Opening ratio
          D(ζ) ~ (Δp > 0) * ((1 - ζ) * Kvo * Δp) + (Δp < 0) * (ζ * Kvc * Δp)
          Aeff_min ~ Mrg * Ann + eps()
          Aeff_max ~ Mst * Ann
          Aeff ~ (Aeff_max - Aeff_min) * ζ + Aeff_min
          # Flow equation
          B ~ (ρ/1000) / (2 * Aeff^2)
          q ~ sqrt(1 / B * abs(Δp)) * sign(Δp)
  end
end

"""
`MynardValve_SemiLunar(; name, ρ, Leff, Mrg, Mst, Ann, Kvc, Kvo)`

Implements the Mynard description for flow across the semilunar valves, full description in [Mynard].
This valve description corresponds to the semilunar valves where interia is an effect we consider.

Note: The minimum level of regurgitation has to be set to machine precision eps()

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)
p is scaled to ensure units are consistent throughout.

Named parameters:

name    name of the element
`ρ`     Blood density in kg/m^3
`Leff`  An effective length in cm
`Mrg`   Level of regurgitation exhibited by a valve in DN
`Mst`   Level of stenosis exhibited by a valve in DN
`Ann`   Annulus area in cm^2
`Kvc`   Valve closing rate coefficent in  cm^2/(dynes*s)
`Kvo`   Valve opening rate coefficent in cm^2/(dynes*s)

p is calculated in mmHg
q is calculated in cm^3/s (ml/s)
"""

@mtkmodel MynardValve_SemiLunar begin
  @extend OnePort()
  @parameters begin
          ρ = ρ_b
          Leff
          Mrg = 0.0
          Mst = 1.0
          Ann
          Kvc
          Kvo
  end
  @variables begin
          Aeff(t)
          ζ(t)
          B(t)
          Aeff_min(t)
          Aeff_max(t)
          L(t)
  end
  begin
          Δp = -mmHg2dynecm2 * Δp # Convert mmHg to dynes/cm^2 = 1 Barye = 1 g/(cm*s^2)
  end
  @equations begin
          # Opening ratio
          D(ζ) ~ (Δp > 0) * ((1 - ζ) * Kvo * Δp) + (Δp < 0) * (ζ * Kvc * Δp)
          Aeff_min ~ Mrg * Ann + eps()
          Aeff_max ~ Mst * Ann
          Aeff ~ (Aeff_max - Aeff_min) * ζ + Aeff_min
          # Flow equation
          B ~ (ρ/1000) / (2 * Aeff^2)
          L ~ (ρ/1000) * Leff / Aeff
          D(q) ~ (Δp - B * q * abs(q)) * 1 / L
  end
end

"""
Junctions
These models represent branches in the arterial network. They are required for the blood gas split.
"""

@mtkmodel Junction2 begin
  @components begin
    in = Pin()
    out1 = Pin()
    out2 = Pin()
  end
  @equations begin
    in.p ~ out1.p
    in.p ~ out2.p
    0 ~ in.q + out1.q + out2.q
    out1.cO₂ ~ instream(in.cO₂)
    out2.cO₂ ~ instream(in.cO₂)
    out1.cCO₂ ~ instream(in.cCO₂)
    out2.cCO₂ ~ instream(in.cCO₂)
  end
end

@mtkmodel Junction3 begin
  @components begin
    in = Pin()
    out1 = Pin()
    out2 = Pin()
    out3 = Pin()
  end
  @equations begin
    in.p ~ out1.p
    in.p ~ out2.p
    in.p ~ out3.p
    0 ~ in.q + out1.q + out2.q + out3.q
    out1.cO₂ ~ instream(in.cO₂)
    out2.cO₂ ~ instream(in.cO₂)
    out3.cO₂ ~ instream(in.cO₂)
    out1.cCO₂ ~ instream(in.cCO₂)
    out2.cCO₂ ~ instream(in.cCO₂)
    out3.cCO₂ ~ instream(in.cCO₂)
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
    RespRate₀ = RespRateₙₒₘ
    IEratio = IE_ratio
    ε = 0.00002         # width of the update pulse
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

"""
Tidal Volume Calculator
This model calculates the tidal volume based on the lung volume and the minimum and maximum lung volumes. It uses a breath hold trigger to update the minimum and maximum lung volumes based on the new detected inspiration/expiration.
"""

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

"""
Lung Gas Exchange
This model represents the gas exchange in the lungs. It is based on the work of Albanese (2016) and gas transfer through the dead space and across the alveoli to the pulmonary capillaries. The long equations significantly slow down the simulation, they are there to avoid the algebraic loop that is created by the D(cppO₂) and D(cppCO₂) terms in the equations:
D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * D(cppO₂))) / V_A
D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * D(cppCO₂))) / V_A
See the standalone SolvingLungGasExchange.jl for a more detailed description of the workflow to remove the algebraic loop..
"""

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
    _DlO₂ = DlO₂ # Diffusion coefficient of O₂ in blood (ml O₂/ml blood/mmHg)
    _DlCO₂ = DlCO₂ # Diffusion coefficient of CO₂ in blood (ml CO₂/ml blood/mmHg)
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

    XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
    XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
    paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
    paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

    SaO₂(t) # O₂ saturation in the arterial blood (%)
  end
  @equations begin
    #### Conservation of mass equations
    D(FDO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FIO₂ - FDO₂), Vr_A * (FDO₂ - FAO₂)) / V_D
    D(FDCO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FICO₂ - FDCO₂), Vr_A * (FDCO₂ - FACO₂)) / V_D

    D(FAO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂) + ((-(ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + ((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))) / (-1 + (-_K*(((_K^2)*((-_p_ws + _pₐₜₘ)^2)*_β₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + (-(_K^2)*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))

    D(FACO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂) - _K*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*(((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-(ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp)) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))) / (-1 + ((_K^2)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*((_K*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-_K*((-_p_ws + _pₐₜₘ)^2)*_β₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))

    Vpp ~ (Vpa - v0pa) * (1 - _sh) + v0pa


    # #### Dissociation equations
    cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
    XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
    cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
    XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

    # #### Instantaneous equilibrium equations
    p_AO₂ ~ pppO₂
    p_ACO₂ ~ pppCO₂

    #### Gas fraction to partial pressure relationships
    p_AO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
    p_ACO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

    #### Mixing between capillary and shunted blood
    caO₂ ~ (1 - _sh) * cppO₂ + _sh * cvO₂
    caCO₂ ~ (1 - _sh) * cppCO₂ + _sh * cvCO₂

    qpp ~ qpa * (1 - _sh)
    qps ~ qpa * _sh

    #### O₂ saturation in arterial blood
    XaO₂ ~ max(((caO₂ / _CₛₐₜO₂) / (1 - (caO₂ / _CₛₐₜO₂))),1e-8)^_h₁
    XaCO₂ ~ max(((caCO₂ / _CₛₐₜCO₂) / (1 - (caCO₂ / _CₛₐₜCO₂))),1e-8)^_h₂

    paO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) + _K₁*XaO₂*_β₂ - _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₂ + 2*_K₂*XaCO₂*_α₂*_β₁)

    paCO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) - _K₁*XaO₂*_β₂ + _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₁ + 2*_K₁*XaO₂*_α₁*_β₂)

    SaO₂ ~ (caO₂ - paO₂ * _sol_O₂) / (_Hgb * _Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
  end
end

"""
Dynamic Heart Power
This model calculates the instantaneous oxygen consumption of the heart based on the left and right ventricle pressures and volumes.
"""

@mtkmodel HeartPower begin
  @parameters begin
    Wh₀ = Whₙₒₘ       # Zero pressure volume (ml)
    MO₂ₙ = MO₂_heart
    τw = τ_heart_w       # Time constant (s)
  end
  @variables begin
    MO₂dyn(t)
    Wh(t)
    wh(t) # Instantaneous power
    Plv(t) # Left ventricle pressure (mmHg)
    Prv(t) # Right ventricle pressure (mmHg)
    DVlv(t) # Left ventricle volume (ml)
    DVrv(t) # Right ventricle volume (ml)
  end
  @equations begin
    MO₂dyn ~ (Wh / Wh₀) * MO₂ₙ
    wh ~ - Plv * DVlv - Prv * DVrv
    D(Wh) ~ (wh - Wh) / τw
  end
end