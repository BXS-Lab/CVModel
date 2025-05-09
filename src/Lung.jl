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
Compartment Definitions
"""

@parameters t

D = Differential(t)

@connector Pin begin
  p(t)
  q(t), [connect = Flow]
end

@mtkmodel DrivenLungPressure begin
  @extend OnePort()
  @parameters begin
    P = p_musmin
    TI = T_I
    TE = T_E
    T0 = Tbreath
    τ = τ_mus
  end
  @variables begin
    t0(t) # Time withing the breathing cycle
    ϕ(t)
  end
  @equations begin
    D(ϕ) ~ 1 / T0
    t0 ~ (ϕ - floor(ϕ)) * T0
    # Δp ~ ifelse(t0 <= TI, 0.01 * t0, 0.1 * t0)
    # Time within the breathing cycle
    Δp ~ ifelse(t0 <= TI, (-P/(TI*TE)*t0^2 + P*T0/(TI*TE)*t0),
          (P/(1-exp(-TE/τ))*(exp(-(t0-TI)/τ)-exp(-TE/τ))))
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
    wall = Pin()
  end
  @parameters begin
    p_ao = 0.0 # mmHg
    R_ml = 1.021 /1000 # mmHg.s/ml
    C_l = 1.27
    V₀_l = 34.4 # ml
    R_lt = 0.3369 /1000 # mmHg.s/ml
    C_tr = 2.38
    V₀_tr = 6.63 # ml
    R_tb = 0.3063 /1000 # mmHg.s/ml
    C_b = 13.1
    V₀_b = 18.7 # ml
    R_bA = 0.0817 /1000 # mmHg.s/ml
    C_A = 200
    # V₀_A = 1.263 # ml
    C_cw = 244.5 # ml/mmHg
    FRC = 2400 # ml
  end
  @variables begin
    P_ao(t)
    P_l(t)
    P_tr(t)
    P_b(t)
    P_A(t)
    P_pl(t)
    V_l(t)
    V_tr(t)
    V_b(t)
    V_A(t)
    V_D(t)
    P_mus(t)
    Vrate_in(t)
    Vrate_A(t)
  end
  begin
    V₀_A = FRC + (C_A * -5.0) - (V₀_l) - (C_tr * 5.0 + V₀_tr) - (C_b * 5.0 + V₀_b)
  end
  @equations begin
    in.p ~ P_ao
    in.q ~ Vrate_in
    wall.q ~ 0
    wall.p ~ P_mus
    C_l * D(P_l) ~ (P_ao - P_l) / R_ml - (P_l - P_tr) / R_lt
    C_tr * (D(P_tr) - D(P_pl)) ~ (P_l - P_tr) / R_lt - (P_tr - P_b) / R_tb
    C_b * (D(P_b) - D(P_pl)) ~ (P_tr - P_b) / R_tb - (P_b - P_A) / R_bA
    C_A * (D(P_A) - D(P_pl)) ~ (P_b - P_A) / R_bA
    C_cw * (D(P_pl) - D(P_mus)) ~ (P_l - P_tr) / R_lt
    Vrate_in ~ (P_ao - P_l) / R_ml
    Vrate_A ~ (P_b - P_A) / R_bA
    V_l ~ C_l * P_l + V₀_l
    V_tr ~ C_tr * (P_tr - P_pl) + V₀_tr
    V_b ~ C_b * (P_b - P_pl) + V₀_b
    V_A ~ C_A * (P_A - P_pl) + V₀_A
    V_D ~ V_l + V_tr + V_b
  end
end

"""
Model Parameters
"""

p_musmin = -5.0 # df* cmH2O2mmHg # Minimum respiratory muscle pressure (mmHg)
RRbreath = 12.0 # Breathing Rate (breaths/min)
IEratio = 0.6 # Inspiratory to Expiratory Ratio
Tbreath = 60.0 / RRbreath # Total Breathing Cycle Time (s)
T_E = Tbreath/(1 + IEratio) # Expiratory Time (s)
T_I = T_E * IEratio # Inspiratory Time (s)
τ_mus = T_E/5 # Respiratory muscle time constant (s)

p₀ = 0.0 # mmHg

"""
Model
"""

@named External = ExternalPressureUB(p_ext=p₀)

@named Breathing = DrivenLungPressure()
@named Lung = Lung()

circ_eqs = [
  connect(External.pext, Lung.in, Breathing.in),
  connect(Lung.wall, Breathing.out)
]

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [
  External, Lung, Breathing # Lung Model Breathing ChestWall
  ])

circ_sys = structural_simplify(circ_model)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_model))

u0 = [
  Lung.P_l => 0.0,
  Lung.P_tr => 0.0,
  Lung.P_b => 0.0,
  Lung.P_A => 0.0,
  Breathing.ϕ => 0.0,
  Lung.P_pl => -5.0
]

prob = ODEProblem(circ_sys, u0, tspan)

@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

p1 = plot(Sol, idxs=[Breathing.out.p], xlims = (0, 15),
        label = "P_mus",
        xlabel = "Time (s)",
        ylabel = "Pressure (cmH2O)")

p2 = plot(Sol, idxs=[Lung.P_pl], xlims = (0, 15),
        label = "P_pl",
        xlabel = "Time (s)",
        ylabel = "Pressure (cmH2O)")

p3 = plot(Sol, idxs=[Lung.P_A], xlims = (0, 15),
        label = "P_A",
        xlabel = "Time (s)",
        ylabel = "Pressure (cmH2O)")

p4 = plot(Sol, idxs=[Lung.Vrate_in], xlims = (0, 15),
        label = "Air Flow",
        xlabel = "Time (s)",
        ylabel = "Flow (ml/s)")

p5 = plot(Sol, idxs=[Lung.V_A + Lung.V_D], xlims = (0, 15),
        label = "V_L",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)")

p6 = plot(Sol, idxs=[Lung.V_A], xlims = (0, 15),
        label = "V_A",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)")

p7 = plot(Sol, idxs=[Lung.V_D], xlims = (0, 15),
        label = "V_D",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)")

display(plot(p1,p2,p3,p4,p5,p6,p7, layout=(4,2), size=(900,600), suptitle="Lungs"))

