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

@mtkmodel CentralChemoreceptors begin
  @parameters begin
    Delay = 1.0
    Gain_cA = 1.0
    Gain_cf = 1.0
    set_point = 1.0
    time_cA = 1.0
    time_cf = 1.0
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
    D(y_A) ~ (-y_A + Gain_cA * (delay.y - set_point)) / time_cA
    D(y_f) ~ (-y_f + Gain_cf * (delay.y - set_point)) / time_cf
  end
end