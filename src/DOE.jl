"""
This file contains the simulation setup and 'Design of Experiments' (DOE) for the cardiovascular model. It currently allows setup of the simulation time, tilt angle protocol, gravity environment, and lower body negative pressure (LBNP) protocol. The @mtkmodel objects can be used to define a range of scenarios for the simulation.
"""

"""
Simulation Time (seconds)
"""

Start_time = 0.0
Time_step = 0.01
Stop_time = 800.0
tspan = (Start_time, Stop_time)

"""
Tilt Angle Protocol (α in radians)
"""

@mtkmodel Alpha begin
  @variables begin
    α(t)
  end

  @parameters begin
    α_min = 0.0 # radians
    t_ramp_start = 200.0   # seconds
    t_ramp_end = 230.0     # seconds
    α_max = π/2        # radians
  end

  @equations begin
    # α ~ 0.0
    α ~ ifelse(t < t_ramp_start, α_min,
         ifelse(t < t_ramp_end,
         α_min + (α_max - α_min) * (t - t_ramp_start) / (t_ramp_end - t_ramp_start),
           α_max))
  end
end

"""
Gravity Environment (g in m/s^2)
"""

@mtkmodel Gravity begin
  @variables begin
    g(t)
  end

  @parameters begin
    constant_gravity = 9.81 # m/s^2
    t_ramp_start = 100.0   # seconds
    t_ramp_end = 130.0     # seconds
    g_min = 1e-6        # m/s^2
  end

  @equations begin
    g ~ constant_gravity
    # g ~ ifelse(t <= t_ramp_start,
    # constant_gravity,
    # ifelse(t >= t_ramp_end,
    #        g_min,
    #        constant_gravity +
    #        (g_min - constant_gravity) * (t - t_ramp_start) / (t_ramp_end - t_ramp_start)))
  end
end

"""
Lower Body Negative Pressure (LBNP) Protocol (p_lbnp in mmHg, should be negative)
"""

@mtkmodel LBNP begin
  @variables begin
    p_lbnp(t)
  end

  @parameters begin
    no_lbnp = 0 # mmHg
    t_ramp_start = 100.0   # seconds
    t_ramp_end = 130.0     # seconds
    lbnp_max = -50        # mmHg
  end

  @equations begin
    p_lbnp ~ no_lbnp
    # p_lbnp ~ ifelse(t <= t_ramp_start,
    # no_lbnp,
    # ifelse(t >= t_ramp_end,
    #        lbnp_max,
    #        no_lbnp +
    #        (lbnp_max - no_lbnp) * (t - t_ramp_start) / (t_ramp_end - t_ramp_start)))
  end
end

"""
Gas Composition Protocol (FIO₂ and FICO₂)
"""

@mtkmodel Atmosphere begin
  @variables begin
    FIO₂(t)
    FICO₂(t)
    pₐₜₘ(t)
    pIO₂(t)
    pICO₂(t)
  end

  @parameters begin
    sl_O₂ = 21.0379/100
    sl_CO₂ = 0.0421/100
    sl_pₐₜₘ = 760 # mmHg
    _p_ws = p_ws # Water vapor pressure (mmHg)
  end

  @equations begin
    FIO₂ ~ sl_O₂
    FICO₂ ~ sl_CO₂
    pₐₜₘ ~ sl_pₐₜₘ
    pIO₂ ~ FIO₂ * (pₐₜₘ - _p_ws) # Partial pressure of O₂ (mmHg)
    pICO₂ ~ FICO₂ * (pₐₜₘ - _p_ws) # Partial pressure of CO₂ (mmHg)
  end
end