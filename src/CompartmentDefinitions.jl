#########################
# Component Elements
#########################

@parameters t

D = Differential(t)

@connector Pin begin
  p(t)
  q(t), [connect = Flow]
end

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

@mtkmodel OnePortWithExtPressure begin
  @components begin
    out = Pin()
    in = Pin()
    ep = Pin()
  end
  @variables begin
    Δp(t)
    q(t)
    pg(t)
  end
  @equations begin
    Δp ~ out.p - in.p
    0 ~ in.q + out.q
    0 ~ ep.q
    q ~ in.q
    pg ~ p - ep.p
  end
end

"""
`Resistor(;name, R=1.0)`

Implements the resistor using Ohm's law to represent a vessels linear resistance to blood flow.

Parameter is in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg
`q` calculated in cm^3/s (ml/s)

Named parameters:

`R`:       Resistance of the vessel to the fluid in mmHg*s/ml
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
`Compliance(; name, V₀=0.0, C=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)`

Implements the compliance of a vessel.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`p` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`V₀`:               Unstressed volume ml

`C`:                Vessel compliance in ml/mmHg


`inP`:             (Bool) formulate in dp/dt (default: false)

`has_ep`:          (Bool) if true, add a parameter `p₀` for pressure offset
                   e.g., for thoracic pressure (default: false)

`p₀`:              External pressure in mmHg (e.g., thorax pressure, default: 0.0)
                   _Note: if this argument is set, it will be used, even if `has_ep` is
                   `false`. `has_ep` only controls if `p₀` will be exposed as a parameter!_

has_variable_ep`: (Bool) expose pin for variable external pressure (default: false)
                   This pin can be connected to another pin or function providing external pressure.
                   _Note: if `has_variable_ep` is set to `true` this pin is created, independent of
                   `has_ep`!_
"""
@component function Compliance(; name, V₀=0.0, C=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)
  @named in = Pin()
  @named out = Pin()

  if has_variable_ep
    @named ep = Pin()
  end

  sts = @variables begin
    V(t)
    p(t)
  end

  ps = @parameters begin
    V₀ = V₀
    C = C
  end

  # Add the thoracic pressure variant

  D = Differential(t)

  eqs = [
    0 ~ in.p - out.p
    p ~ in.p
  ]

  if has_variable_ep
    push!(sts,
      (@variables p_rel(t))[1]
    )
    if has_ep
      push!(ps,
        (@parameters p₀ = p₀)[1]
      )
    end
    push!(eqs,
      p_rel ~ ep.p + p₀,
      ep.q ~ 0
    )
  elseif has_ep
    push!(ps,
      (@parameters p₀ = p₀)[1]
    )
    p_rel = p₀
  else
    p_rel = p₀
  end

  if inP
    push!(eqs,
      V ~ (p - p_rel) * C + V₀,
      D(p) ~ (in.q + out.q) * 1 / C
    )
  else
    push!(eqs,
      p ~ (V - V₀) / C + p_rel,
      D(V) ~ in.q + out.q
    )
  end

  if has_variable_ep
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out, ep)
  else
    compose(ODESystem(eqs, t, sts, ps; name=name), in, out)
  end
end

"""
`HeldtChamber(;name, V₀, p₀=0.0, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift=0.0)`

Implemention of a ventricle following Heldt.

This model uses external helper function `HeldtElastance`
which describes the elastance function.

Named parameters:

`V₀`     stress-free volume (zero pressure volume)

`p₀`     pressure offset (defaults to zero)
         this is present in the original paper, so is
         provided here for conformity. Defaults to 0.0

`Eₘᵢₙ`   minimum elastance

`τ`      pulse length

`τₑₛ`    end systolic time (end of rising cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle

`inP`:    (Bool) formulate in dp/dt (default: false)
"""
@mtkmodel HeldtChamber begin
  @structural_parameters begin
    inP = false
  end
  @variables begin
    V(t)
    p(t)
  end
  @parameters begin
    V₀
    p₀ = 0.0
    Eₘᵢₙ
    Eₘₐₓ
    τ
    τₑₛ
    Eshift
  end

  @components begin
    in = Pin()
    out = Pin()
  end
  begin
    E = HeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)
    DE = DHeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)
    p_rel = p₀
  end

  @equations begin
    0 ~ in.p - out.p
    p ~ in.p
    if inP

      V ~ (p - p_rel) / E + V₀
      D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE
    else
      p ~ (V - V₀) * E + p_rel
      D(V) ~ in.q + out.q
    end
  end
end


"""
`HeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)`

Elastance function `E(t)` for ventricle simulation based on Heldt's
double cosine function.

Parameters:

`Eₘᵢₙ`: minimum elastance (diastole)

`Eₘₐₓ`: maximum elastance (systole)

`τₑₛ`: end systolic time (end of rising cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle
"""
function HeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)

  tᵢ = rem(t + (1 - Eshift) * τ, τ)

  Eₚ = (tᵢ <= τₑₛ) * (1 - cos(tᵢ / τₑₛ * pi)) / 2 +
       (tᵢ > τₑₛ) * (tᵢ <= (1.5 * τₑₛ)) * (1 + cos((tᵢ - τₑₛ) / (τₑₛ) * 2 * pi)) / 2 +
       (tᵢ <= (1.5 * τₑₛ)) * 0

  E = Eₘᵢₙ + (Eₘₐₓ - Eₘᵢₙ) * Eₚ

  return E
end


"""
DHeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)

Helper function for `ShiChamber`

Derivative of the elastance function `E(t)` for ventricle simulation based on Shi's
double cosine function.

Parameters:

`Eₘᵢₙ`: minimum elastance (diastole)

`Eₘₐₓ`: maximum elastance (systole)

`τₑₛ`: end systolic time (end of rising cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle
"""
function DHeldtElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, Eshift)

  tᵢ = rem(t + (1 - Eshift) * τ, τ)

  DEₚ = (tᵢ <= τₑₛ) * pi / τₑₛ * sin(tᵢ / τₑₛ * pi) / 2 +
        (tᵢ > τₑₛ) * (tᵢ <= (1.5 * τₑₛ)) * pi / τₑₛ * -sin((τₑₛ - tᵢ) / (τₑₛ) * 2 * pi) +
        (tᵢ <= (1.5 * τₑₛ)) * 0
  DE = (Eₘₐₓ - Eₘᵢₙ) * DEₚ

  return DE
end

"""
`ResistorDiode(;name, R=1e-3)`

Implements the resistance across a valve following Ohm's law exhibiting diode like behaviour.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)

Named parameters:

`R`     Resistance across the valve in mmHg*s/ml
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
`Inductance(;name, L=1.0)`

Implements the inductance to represent blood inertance.

Parameters are in the cm, g, s system.
Pressures in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`L`:       Inertia of the fluid in mmHg*s^2/ml
"""
@mtkmodel Inductance begin
  @extend OnePort()
  @parameters begin
    L = 1.0
  end
  @equations begin
    D(q) ~ -Δp / L
  end
end

#########################
# Define Lumped Compartments
#########################

"""
`Systemic(;name, R=1.0, C=1.0)`

Implements the compliance, resistor systemic compartment.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`R`:       Component resistance in mmHg*s/ml

`C`:       Component compliance in ml/mmHg

"""
@mtkmodel StandardComp begin
  @structural_parameters begin
    R = 1.0
    C = 1.0
    p₀=0.0
    V₀=1.0
  end
  @variables begin
    Δp(t)
    q(t)
  end
  @components begin
    in = Pin()
    out = Pin()
    R = Resistor(R=R)
    C = Compliance(V₀=V₀, C=C, inP=false, has_ep=true, has_variable_ep=false, p₀=p₀)
  end
  @equations begin
    Δp ~ out.p - in.p
    q ~ in.q
    connect(in, C.in)
    connect(C.out, R.in)
    connect(R.out, out)
  end
end

"""
`CRL(;name, C=1.0, R=1.0, L=1.0)`

Implements the compliance, resistor, inductance subsystem.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`C`:       Component compliance in ml/mmHg

`R`:       Component resistance in mmHg*s/ml

`L`:       Component blood inertia in mmHg*s^2/ml
"""
@mtkmodel CRL begin
  @structural_parameters begin
    C = 1.0
    R = 1.0
    L = 1.0
    p₀= 0.0
    V₀= 1.0
  end
  @variables begin
    Δp(t)
    q(t)
  end
  @components begin
    in = Pin()
    out = Pin()
    C = Compliance(V₀=V₀, C=C, inP=false, has_ep=true, has_variable_ep=false, p₀=p₀)
    R = Resistor(R=R)
    L = Inductance(L=L)
  end
  @equations begin
    Δp ~ out.p - in.p
    q ~ in.q
    connect(in, C.in)
    connect(C.out, R.in)
    connect(R.out, L.in)
    connect(L.out, out)
  end
end