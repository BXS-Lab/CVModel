"""
STANDALONE FILE – NOT CONNECTED TO THE MAIN CVModel.jl
This file is a standalone simulation of the impulse response of each branch of the reflex control system. You can use this to test the impulse response and vary the order of the Pade approximation. The code is structurally identical to the main CVModel.jl file, but it is not connected to the rest of the model.
"""

"""
Preamble
This code section sets up the environment for this standalone simulation and defines the necessary model parameters. In the main CVModel.jl, these are all defined elsewhere. Here you can vary the order of the Pade approximation to see the effect on the computational time.
"""

using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using Symbolics

Start_time = 0.0
Time_step = 0.01
Stop_time = 200
tspan = (Start_time, Stop_time)

reflex_delay_order = 10 # Order of the Pade delay
reflex_delay_init = zeros(reflex_delay_order) # Initial condition for the delay

abr_αr_delay = 2.5
abr_αr_peak = 3.5
abr_αr_end = 30.0

abr_αv_delay = 5.0
abr_αv_peak = 10.0
abr_αv_end = 42.0

abr_β_delay = 2.5
abr_β_peak = 3.5
abr_β_end = 15.0

abr_para_delay = 0.59
abr_para_peak = 0.70
abr_para_end = 1.0

cpr_αr_delay = 2.5
cpr_αr_peak = 5.5
cpr_αr_end = 35.0

cpr_αv_delay = 5.0
cpr_αv_peak = 9.0
cpr_αv_end = 40.0

@parameters t
D = Differential(t)

"""
Define Compartments
This section defines the compartments for the transfer function. Their implementation is identical to the Reflex.jl file.
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
Driver for the Impulse Response
This code section creates a simple driver that generates a dirac delta function to the system in order to output the impulse response.
"""

@mtkmodel ImpulseDriver begin
  @parameters begin
    magnitude = 100.0  # Height of the impulse
    width = 0.01       # Duration of the impulse
  end
  @variables begin
    u(t)
  end
  @equations begin
    u ~ ifelse(t < width, magnitude, 0.0)
  end
end

"""
Instance and Connect the Compartments
Here we instance the compartments for this standalone simulation and connect them all to the impulse driver in parallel. Each instance represents a different part of the reflex control: Arterial Baroreflex (α-Sympathetic to Resistance, α-Sympathetic to Volume, β-Sympathetic, and Parasympathetic) and Cardiopulmonary Reflex (α-Sympathetic to Resistance and α-Sympathetic to Volume).
"""

@named abr_αr = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_αr_delay, reflex_peak = abr_αr_peak, reflex_end = abr_αr_end)
@named abr_αv = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_αv_delay, reflex_peak = abr_αv_peak, reflex_end = abr_αv_end)
@named abr_β = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_β_delay, reflex_peak = abr_β_peak, reflex_end = abr_β_end)
@named abr_para = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_para_delay, reflex_peak = abr_para_peak, reflex_end = abr_para_end)
@named cpr_αr = TransferFunction(delay_order = reflex_delay_order, reflex_delay = cpr_αr_delay, reflex_peak = cpr_αr_peak, reflex_end = cpr_αr_end)
@named cpr_αv = TransferFunction(delay_order = reflex_delay_order, reflex_delay = cpr_αv_delay, reflex_peak = cpr_αv_peak, reflex_end = cpr_αv_end)
@named input = ImpulseDriver()

circ_eqs = [
  input.u ~ abr_αr.u,
  input.u ~ abr_αv.u,
  input.u ~ abr_β.u,
  input.u ~ abr_para.u,
  input.u ~ cpr_αr.u,
  input.u ~ cpr_αv.u
]

"""
Compose and Solve
Here we compose the system, initialize it to zero throughout, and solve it. The code is structurally identical to that found in the main CVModel.jl file.
"""

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [input, abr_αr, abr_αv, abr_β, abr_para, cpr_αr, cpr_αv])

circ_sys = structural_simplify(circ_model)

equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_sys))

u0 = [
  abr_αr.tfdelay.tftime.x => reflex_delay_init,
  abr_αr.tfdelay.double_integrator.v => 0.0,
  abr_αr.tfdelay.double_integrator.y => 0.0,
  abr_αr.tfpeak.tftime.x => reflex_delay_init,
  abr_αr.tfpeak.double_integrator.v => 0.0,
  abr_αr.tfpeak.double_integrator.y => 0.0,
  abr_αr.tfend.tftime.x => reflex_delay_init,
  abr_αr.tfend.double_integrator.v => 0.0,
  abr_αr.tfend.double_integrator.y => 0.0,
  abr_αv.tfdelay.tftime.x => reflex_delay_init,
  abr_αv.tfdelay.double_integrator.v => 0.0,
  abr_αv.tfdelay.double_integrator.y => 0.0,
  abr_αv.tfpeak.tftime.x => reflex_delay_init,
  abr_αv.tfpeak.double_integrator.v => 0.0,
  abr_αv.tfpeak.double_integrator.y => 0.0,
  abr_αv.tfend.tftime.x => reflex_delay_init,
  abr_αv.tfend.double_integrator.v => 0.0,
  abr_αv.tfend.double_integrator.y => 0.0,
  abr_β.tfdelay.tftime.x => reflex_delay_init,
  abr_β.tfdelay.double_integrator.v => 0.0,
  abr_β.tfdelay.double_integrator.y => 0.0,
  abr_β.tfpeak.tftime.x => reflex_delay_init,
  abr_β.tfpeak.double_integrator.v => 0.0,
  abr_β.tfpeak.double_integrator.y => 0.0,
  abr_β.tfend.tftime.x => reflex_delay_init,
  abr_β.tfend.double_integrator.v => 0.0,
  abr_β.tfend.double_integrator.y => 0.0,
  abr_para.tfdelay.tftime.x => reflex_delay_init,
  abr_para.tfdelay.double_integrator.v => 0.0,
  abr_para.tfdelay.double_integrator.y => 0.0,
  abr_para.tfpeak.tftime.x => reflex_delay_init,
  abr_para.tfpeak.double_integrator.v => 0.0,
  abr_para.tfpeak.double_integrator.y => 0.0,
  abr_para.tfend.tftime.x => reflex_delay_init,
  abr_para.tfend.double_integrator.v => 0.0,
  abr_para.tfend.double_integrator.y => 0.0,
  cpr_αr.tfdelay.tftime.x => reflex_delay_init,
  cpr_αr.tfdelay.double_integrator.v => 0.0,
  cpr_αr.tfdelay.double_integrator.y => 0.0,
  cpr_αr.tfpeak.tftime.x => reflex_delay_init,
  cpr_αr.tfpeak.double_integrator.v => 0.0,
  cpr_αr.tfpeak.double_integrator.y => 0.0,
  cpr_αr.tfend.tftime.x => reflex_delay_init,
  cpr_αr.tfend.double_integrator.v => 0.0,
  cpr_αr.tfend.double_integrator.y => 0.0,
  cpr_αv.tfdelay.tftime.x => reflex_delay_init,
  cpr_αv.tfdelay.double_integrator.v => 0.0,
  cpr_αv.tfdelay.double_integrator.y => 0.0,
  cpr_αv.tfpeak.tftime.x => reflex_delay_init,
  cpr_αv.tfpeak.double_integrator.v => 0.0,
  cpr_αv.tfpeak.double_integrator.y => 0.0,
  cpr_αv.tfend.tftime.x => reflex_delay_init,
  cpr_αv.tfend.double_integrator.v => 0.0,
  cpr_αv.tfend.double_integrator.y => 0.0
]

prob = ODEProblem(circ_sys, u0, tspan)
@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

"""
Plot
This line plots the impulse response for each reflex type.
"""

display(plot(Sol, idxs=[abr_αr.y,abr_αv.y,abr_β.y,abr_para.y,cpr_αr.y,cpr_αv.y], xlims=(0, 60), ylims=(0, 0.2)))

"""
Debug Script – Unit Impulse
This code section extracts the impulse response for each reflex and numerically integrates it to find the area under the curve. The areas are printed to the console. By increasing the order of the Pade approximation, you can improve the accuracy of the unit impulse at the expense of computational time.
"""

ts = Sol.t
ys_abr_αr = Sol[abr_αr.y, :]
ys_abr_αv = Sol[abr_αv.y, :]
ys_abr_β = Sol[abr_β.y, :]
ys_abr_para = Sol[abr_para.y, :]
ys_cpr_αr = Sol[cpr_αr.y, :]
ys_cpr_αv = Sol[cpr_αv.y, :]

# Calculate area for each response
area_abr_αr = sum(diff(ts) .* (ys_abr_αr[1:end-1] .+ ys_abr_αr[2:end]) ./ 2)
area_abr_αv = sum(diff(ts) .* (ys_abr_αv[1:end-1] .+ ys_abr_αv[2:end]) ./ 2)
area_abr_β = sum(diff(ts) .* (ys_abr_β[1:end-1] .+ ys_abr_β[2:end]) ./ 2)
area_abr_para = sum(diff(ts) .* (ys_abr_para[1:end-1] .+ ys_abr_para[2:end]) ./ 2)
area_cpr_αr = sum(diff(ts) .* (ys_cpr_αr[1:end-1] .+ ys_cpr_αr[2:end]) ./ 2)
area_cpr_αv = sum(diff(ts) .* (ys_cpr_αv[1:end-1] .+ ys_cpr_αv[2:end]) ./ 2)

println("Areas under curves:")
println("ABR αr: ", area_abr_αr)
println("ABR αv: ", area_abr_αv)
println("ABR β: ", area_abr_β)
println("ABR para: ", area_abr_para)
println("CPR αr: ", area_cpr_αr)
println("CPR αv: ", area_cpr_αv)