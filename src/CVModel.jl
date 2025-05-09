"""
Cardiovascular Model
Version 3.0 (May 3rd, 2025)
BXS Lab, UC Davis; Authors: RS Whittle, AJ Kondoor, HS Vellore
Contact info: Dr. Rich Whittle – Department of Mechanical and Aerospace Engineering, UC Davis, One Shields Ave, Davis CA 95616 (rswhittle@ucdavis.edu)
This model is a simulation of the human cardiovascular system, including a four chamber heart, arteries, veins, and microcirculation. It includes reflex control for the arterial baroreflex (ABR) and cardiopulmonary reflex (CPR), as well as hydrostatic effects, interstitial fluid flow, and external tissue pressures. The model is designed to simulate various physiological scenarios, including a tilt angle protocol, altered-gravity environment, and lower body negative pressure (LBNP) protocol. The underlying equations are based on the work of Heldt (2004), Zamanian (2007), Diaz Artiles (2015), and Whittle (2023). The model is implemented using the ModelingToolkit.jl package in Julia.
"""
### TODO: (1) Collapse in vessels (2) Bring Intrathoracic Pressure into Lung Model

module CVModel
display("Cardiovascular Model v3.0 (May 3rd, 2025) - BXS Lab, UC Davis")

"""
Preamble
These are the required Julia packages.
"""

using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using Symbolics

"""
Model Parameters
This file loads all of the model parameters.
"""

include("ModelParameters.jl")
using .ModelParams

"""
Compartment Definitions
These two files contain the definitions of the compartments used in the model. "CompartmentDefinitions.jl" contains the definitions of the hemodynamic compartments, while "reflex.jl" contains the definitions of the reflex arcs (ABR and CPR).
"""

include("CompartmentDefinitions.jl")
include("Reflex.jl")

"""
Design of Experiments
This file contains the simulation setup and 'Design of Experiments' (DOE). It currently allows setup of the simulation time, tilt angle protocol, gravity environment, and lower body negative pressure (LBNP) protocol.
"""

include("DOE.jl")

"""
Find Initial Conditions
This file uses Heldt's method to find the initial pressure vector for the model using a modified Newton-Raphson algorithm. It also calculates the initial interstitial volume (modeled as an equilibrium state) based on the starting angle, gravity, and LBNP.
"""

include("Initial.jl")

"""
Instance Compartments
This section of code instances the compartments used in the model, based on the definitions created above and the model parameters.
"""

#### Sino-Atrial Node
@named SA = SANode(RR₀=RRₙₒₘ, τₐᵥ = τₐᵥ)

#### Right Heart
@named RA = HeldtChamber(V₀=v0_ra, Eₘᵢₙ=Ed_ra, Eₘₐₓ=Ees_ra, τₑₛ=τₐₛ, inP=true)
@named R_tricuspid = MynardValve_Atrioventricular(Ann=Ann_tv, Kvc=Kvc_tv, Kvo=Kvo_tv)
@named RV = HeldtChamber(V₀=v0_rv, Eₘᵢₙ=Ed_rv, Eₘₐₓ=Ees_rv, Elimit = Elimit_rv, τₑₛ=τᵥₛ, inP=true, has_abr=true)
@named R_pulmonary = MynardValve_SemiLunar(Leff=Leff_pv, Ann=Ann_pv, Kvc=Kvc_pv, Kvo=Kvo_pv)

#### Pulmonary Circulation
@named Pulm_art = Artery(C=Cpa, R=Rpa, V₀=v0pa, has_hydrostatic=false, has_tissue=false, has_inertia=false)
@named Pulm_cap = StarlingResistor(R=Rpc, h=h_Lungs)
@named Pulm_vein = Vein(C=Cpv, R=Rpv, V₀=v0pv, has_hydrostatic=false, has_tissue=false)

#### Left Heart
@named LA = HeldtChamber(V₀=v0_la, Eₘᵢₙ=Ed_la, Eₘₐₓ=Ees_la, τₑₛ=τₐₛ, inP=true)
@named R_mitral = MynardValve_Atrioventricular(Ann=Ann_mv, Kvc=Kvc_mv, Kvo=Kvo_mv)
@named LV = HeldtChamber(V₀=v0_lv, Eₘᵢₙ=Ed_lv, Eₘₐₓ=Ees_lv, Elimit = Elimit_lv, τₑₛ=τᵥₛ, inP=true, has_abr=true)
@named R_aortic = MynardValve_SemiLunar(Leff=Leff_av, Ann=Ann_av, Kvc=Kvc_av, Kvo=Kvo_av)

#### Coronary Circulation
@named Cor_art = Artery(C=Cca, R=Rca, V₀=v0ca, has_hydrostatic=false, has_tissue=false)
@named Cor_cap = VarResistor()
@named Cor_vein = Vein(C=Ccv, R=Rcv, V₀=v0cv, has_hydrostatic=false, has_tissue=false)

#### Arterial Circulation
@named Asc_A = Artery(C=C_Asc_A, R=R_Asc_A, V₀=v0_Asc_A, h=h_Asc_A, has_tissue=false, has_inertia=false)
@named BC_A = Artery(C=C_BC_A, R=R_BC_A, V₀=v0_BC_A, h=h_BC_A, has_tissue=false, L=L_BC_A)
@named UpBd_art = Artery(C=C_UpBd_art, R=R_UpBd_art, V₀=v0_UpBd_art, h=h_UpBd_art, rad=rad_UB, L=L_UpBd_art)
@named Thor_A = Artery(C=C_Thor_A, R=R_Thor_A, V₀=v0_Thor_A, h=h_Thor_A, has_tissue=false, L=L_Thor_A)
@named Abd_A = Artery(C=C_Abd_A, R=R_Abd_A, V₀=v0_Abd_A, h=h_Abd_A, rad=rad_Abd, L=L_Abd_A)
@named Renal_art = Artery(C=C_Renal_art, R=R_Renal_art, V₀=v0_Renal_art, h=h_Renal_art, rad=rad_Abd, L=L_Renal_art)
@named Splanchnic_art = Artery(C=C_Splanchnic_art, R=R_Splanchnic_art, V₀=v0_Splanchnic_art, h=h_Splanchnic_art, rad=rad_Abd, L=L_Splanchnic_art)
@named Leg_art = Artery(C=C_Leg_art, R=R_Leg_art, V₀=v0_Leg_art, h=h_Leg_art, con=con_Leg_art, rad=rad_Leg, L=L_Leg_art)

#### Venous Circulation
@named UpBd_vein = Vein(C=C_UpBd_vein, R=R_UpBd_vein, V₀=v0_UpBd_vein, h=h_UpBd_vein, has_valve=true, has_abr=true, has_cpr=true, rad=rad_UB)
@named SVC = Vein(C=C_SVC, R=R_SVC, V₀=v0_SVC, h=h_SVC, has_tissue=false)
@named Renal_vein = Vein(C=C_Renal_vein, R=R_Renal_vein, V₀=v0_Renal_vein, h=h_Renal_vein, has_abr=true, has_cpr=true, rad=rad_Abd, V_min=vₘᵢₙ_Renal_vein)
@named Splanchnic_vein = Vein(C=C_Splanchnic_vein, R=R_Splanchnic_vein, V₀=v0_Splanchnic_vein, h=h_Splanchnic_vein, is_nonlinear=true, Flow_div = Flow_Splanchnic_vein, V_max=vM_Splanchnic_vein, has_abr=true, has_cpr=true, rad=rad_Abd)
@named Leg_vein = Vein(C=C_Leg_vein, R=R_Leg_vein, V₀=v0_Leg_vein, h=h_Leg_vein, con=con_Leg_vein, has_valve=true, is_nonlinear=true, Flow_div = Flow_Leg_vein, V_max=vM_Leg_vein, has_abr=true, has_cpr=true, rad=rad_Leg)
@named Abd_veins = Vein(C=C_Abd_veins, R=R_Abd_veins, V₀=v0_Abd_veins, h=h_Abd_veins, is_nonlinear=true, Flow_div = Flow_Abd_veins, V_max=vM_Abd_vein, rad=rad_Abd)
@named Thor_IVC = Vein(C=C_Thor_IVC, R=R_Thor_IVC, V₀=v0_Thor_IVC, h=h_Thor_IVC, has_tissue=false)

#### Microcirculation
@named UpBd_cap = VarResistor()
@named Renal_cap = VarResistor()
@named Splanchnic_cap = VarResistor()
@named Leg_cap = VarResistor()

#### Interstitial Compartment
@named Interstitial = InterstitialCompartment(Vmtilt=Vmax_tilt, Vmlbnp=Vmax_lbnp, τ=τint)

#### External Pressures
@named Intrathoracic = IntrathoracicPressure()
@named Abdominal = IntraAbdominalPressure(p_abd=p_abd)
@named External = ExternalPressureUB(p_ext=p₀)
@named ExternalLBNP = ExternalPressureLB(p_ext=p₀)

#### Design of Experiments
# These are the drivers to implement the protocols defined in the "doe.jl" file.
@named alpha_driver = Alpha()
@named gravity_driver = Gravity()
@named lbnp_driver = LBNP()

#### Reflex Arcs
# There are two afferents, one for each arc (ABR and CPR). The separate effects (e.g., α-Sympathetic, β-Sympathetic, Parasympathetic) are modeled as separate transfer functions.
@named CPRafferent = Afferent(p_set=p_cpr, gain=gain_cpr)
@named ABRafferent = Afferent(p_set=p_abr, gain=gain_abr)

@named abr_αr = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_αr_delay, reflex_peak = abr_αr_peak, reflex_end = abr_αr_end)
@named abr_αv = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_αv_delay, reflex_peak = abr_αv_peak, reflex_end = abr_αv_end)
@named abr_β = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_β_delay, reflex_peak = abr_β_peak, reflex_end = abr_β_end)
@named abr_para = TransferFunction(delay_order = reflex_delay_order, reflex_delay = abr_para_delay, reflex_peak = abr_para_peak, reflex_end = abr_para_end)
@named cpr_αr = TransferFunction(delay_order = reflex_delay_order, reflex_delay = cpr_αr_delay, reflex_peak = cpr_αr_peak, reflex_end = cpr_αr_end)
@named cpr_αv = TransferFunction(delay_order = reflex_delay_order, reflex_delay = cpr_αv_delay, reflex_peak = cpr_αv_peak, reflex_end = cpr_αv_end)

#### Lung model
@named RespMuscles = RespiratoryMuscles()
@named Lung = Lung()

"""
Structural Connections
This section of code connects the instanced compartments together to form the cardiovascular system. The hemodynamic connections are made using the "connect" function, which connects pins linked in pressure and flow via Kirchhoff's laws. The DOE integrations and reflex connections are just included as direct signal connections to the respective compartments.
"""

circ_eqs = [
  #### Heart and Pulmonary System
  connect(RA.out, R_tricuspid.in),
  connect(R_tricuspid.out, RV.in),
  connect(RV.out, R_pulmonary.in),
  connect(R_pulmonary.out,Pulm_art.in),
  connect(Pulm_art.out, Pulm_cap.in),
  connect(Pulm_cap.out, Pulm_vein.in),
  connect(Pulm_vein.out, LA.in),
  connect(LA.out, R_mitral.in),
  connect(R_mitral.out, LV.in),
  connect(LV.out, R_aortic.in),
  connect(R_aortic.out, Asc_A.in),

  #### Arterial Tree
  connect(Asc_A.out, BC_A.in, Thor_A.in, Cor_art.in),
  connect(BC_A.out, UpBd_art.in),
  connect(Thor_A.out, Abd_A.in),
  connect(Abd_A.out, Renal_art.in, Splanchnic_art.in, Leg_art.in),

  #### Coronary Circulation
  connect(Cor_art.out, Cor_cap.in),
  connect(Cor_cap.out, Cor_vein.in),

  #### Upper Body Circulation
  connect(UpBd_art.out, UpBd_cap.in),
  connect(UpBd_cap.out, UpBd_vein.in),

  #### Renal Circulation
  connect(Renal_art.out, Renal_cap.in),
  connect(Renal_cap.out, Renal_vein.in),

  #### Splanchnic Circulation
  connect(Splanchnic_art.out, Splanchnic_cap.in),
  connect(Splanchnic_cap.out, Splanchnic_vein.in),

  #### Leg Circulation
  connect(Leg_art.out, Leg_cap.in),
  connect(Leg_cap.out, Leg_vein.in),

  #### Venous Tree
  connect(UpBd_vein.out, SVC.in),
  connect(Abd_veins.in, Renal_vein.out, Splanchnic_vein.out, Leg_vein.out),
  connect(Abd_veins.out, Thor_IVC.in),
  connect(Thor_IVC.out, SVC.out, Cor_vein.out, RA.in),

  #### External Pressures
  connect(Intrathoracic.pth, Asc_A.ep, BC_A.ep, Thor_A.ep, SVC.ep, Thor_IVC.ep, RA.ep, RV.ep, Pulm_art.ep, Pulm_vein.ep, LA.ep, LV.ep, Cor_art.ep, Cor_vein.ep),
  connect(Abdominal.pabd, Abd_A.ep, Renal_art.ep, Splanchnic_art.ep, Renal_vein.ep, Splanchnic_vein.ep, Abd_veins.ep),
  connect(External.pext, UpBd_art.ep, UpBd_vein.ep, Lung.in, RespMuscles.in),
  connect(ExternalLBNP.pext, Leg_art.ep, Leg_vein.ep),

  connect(Lung.chestwall, RespMuscles.out),

  Pulm_cap.pₐₗᵥ ~ Lung.p_A,
  Intrathoracic.pₚₗ ~ Lung.pₚₗ,

  #### Interstitial Connections (Direct Connections)
  Splanchnic_vein.C.qint ~ Interstitial.Qint,
  Leg_vein.C.qint ~ Interstitial.Qint,
  Abd_veins.C.qint ~ Interstitial.Qint,

  #### Tilt Equations (Direct Connections)
  Asc_A.α ~ alpha_driver.α,
  BC_A.α ~ alpha_driver.α,
  UpBd_art.α ~ alpha_driver.α,
  Thor_A.α ~ alpha_driver.α,
  Abd_A.α ~ alpha_driver.α,
  Renal_art.α ~ alpha_driver.α,
  Splanchnic_art.α ~ alpha_driver.α,
  Leg_art.α ~ alpha_driver.α,
  UpBd_vein.α ~ alpha_driver.α,
  SVC.α ~ alpha_driver.α,
  Renal_vein.α ~ alpha_driver.α,
  Splanchnic_vein.α ~ alpha_driver.α,
  Leg_vein.α ~ alpha_driver.α,
  Abd_veins.α ~ alpha_driver.α,
  Thor_IVC.α ~ alpha_driver.α,
  Interstitial.α ~ alpha_driver.α,
  Intrathoracic.α ~ alpha_driver.α,
  Pulm_art.α ~ alpha_driver.α,
  Pulm_vein.α ~ alpha_driver.α,
  Pulm_cap.α ~ alpha_driver.α,
  Lung.α ~ alpha_driver.α,

  #### Gravity Equations (Direct Connections)
  Asc_A.g ~ gravity_driver.g,
  BC_A.g ~ gravity_driver.g,
  UpBd_art.g ~ gravity_driver.g,
  Thor_A.g ~ gravity_driver.g,
  Abd_A.g ~ gravity_driver.g,
  Renal_art.g ~ gravity_driver.g,
  Splanchnic_art.g ~ gravity_driver.g,
  Leg_art.g ~ gravity_driver.g,
  UpBd_vein.g ~ gravity_driver.g,
  SVC.g ~ gravity_driver.g,
  Renal_vein.g ~ gravity_driver.g,
  Splanchnic_vein.g ~ gravity_driver.g,
  Leg_vein.g ~ gravity_driver.g,
  Abd_veins.g ~ gravity_driver.g,
  Thor_IVC.g ~ gravity_driver.g,
  Interstitial.g ~ gravity_driver.g,
  Intrathoracic.g ~ gravity_driver.g,
  Pulm_art.g ~ gravity_driver.g,
  Pulm_vein.g ~ gravity_driver.g,
  Pulm_cap.g ~ gravity_driver.g,
  Lung.g ~ gravity_driver.g,

  #### LBNP Equations (Direct Connections)
  ExternalLBNP.p_lbnp ~ lbnp_driver.p_lbnp,
  Interstitial.p_lbnp ~ lbnp_driver.p_lbnp,

  #### Reflex Arc Afferent Inputs (Sensed Pressures)
  CPRafferent.e ~ (RA.pₜₘ),
  ABRafferent.e ~ (Asc_A.C.pₜₘ+(BC_A.C.pₜₘ + (ρ_b*gravity_driver.g*(h_cs/100)*sin(alpha_driver.α)*Pa2mmHg)))/2.0,

  #### Reflex Arc Afferent Outputs Connected to Transfer Function Inputs
  ABRafferent.δ ~ abr_αr.u,
  ABRafferent.δ ~ abr_αv.u,
  ABRafferent.δ ~ abr_β.u,
  ABRafferent.δ ~ abr_para.u,
  CPRafferent.δ ~ cpr_αr.u,
  CPRafferent.δ ~ cpr_αv.u,

  #### Reflex Arc Efferents
  # The outputs of the transfer functions are connected directly to the relevant compartments via static gains defined in the parameters file.

  # Coronary Resistance
  Cor_cap.R ~ Rcc,

  # Arterial Resistance (ABR and CPR)
  UpBd_cap.R ~ R_UpBd_cap + (Gabr_r * abr_αr.y) + (Gcpr_r * cpr_αr.y),
  Renal_cap.R ~ R_Renal_cap + (Gabr_r * abr_αr.y) + (Gcpr_r * cpr_αr.y),
  Splanchnic_cap.R ~ R_Splanchnic_cap + (Gabr_r * abr_αr.y) + (Gcpr_r * cpr_αr.y),
  Leg_cap.R ~ R_Leg_cap + (Gabr_r * abr_αr.y) + (Gcpr_r * cpr_αr.y),

  # Venous Tone (ABR and CPR)
  UpBd_vein.C.Vabr ~ Gabr_vub * abr_αv.y,
  Renal_vein.C.Vabr ~ Gabr_vre * abr_αv.y,
  Splanchnic_vein.C.Vabr ~ Gabr_vsp * abr_αv.y,
  Leg_vein.C.Vabr~ Gabr_vlb * abr_αv.y,
  UpBd_vein.C.Vcpr ~ Gcpr_vub * cpr_αv.y,
  Renal_vein.C.Vcpr ~ Gcpr_vre * cpr_αv.y,
  Splanchnic_vein.C.Vcpr ~ Gcpr_vsp * cpr_αv.y,
  Leg_vein.C.Vcpr ~ Gcpr_vlb * cpr_αv.y,

  # Ventricular Contractility (ABR)
  SA.Eabr_rv ~ Gabr_erv * abr_β.y,
  SA.Eabr_lv ~ Gabr_elv * abr_β.y,

  # Heart Rate (ABR)
  SA.RRabr ~ (Gabr_rrsymp * abr_β.y) + (Gabr_rrpara * abr_para.y),

  #### Sino-Atrial Node Connections
  # The SA node holds the RR Interval and Contractility adjustments, only updating at the start of a new cardiac cycle. The φ signal is the modulated cardiac cycle (atria are offset to contract before the ventricles). The τ signal is the beat-held instantaneous RR interval.
  RV.Eabr_held ~ SA.Eabr_rv_held,
  LV.Eabr_held ~ SA.Eabr_lv_held,
  RA.ϕ ~ SA.ϕ_wrapped_atria,
  RV.ϕ ~ SA.ϕ_wrapped,
  LA.ϕ ~ SA.ϕ_wrapped_atria,
  LV.ϕ ~ SA.ϕ_wrapped,
  RA.τ ~ SA.RR_held,
  RV.τ ~ SA.RR_held,
  LA.τ ~ SA.RR_held,
  LV.τ ~ SA.RR_held
]

"""
Define ODEs
This section of the code composes the system of ordinary differential equations (ODEs) using the previously defined compartments and equations. The structural simplify is key to reducing the complexity of the system. Note: a common error is to not include all instanced compartments in the compose function.
"""

@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [SA, RA, R_tricuspid, RV, R_pulmonary, Pulm_art, Pulm_cap, Pulm_vein, LA, R_mitral, LV, R_aortic, # Heart and Lungs
  Cor_art, Cor_cap, Cor_vein, # Coronary Circulation
  Asc_A, BC_A, UpBd_art, Thor_A, Abd_A, Renal_art, Splanchnic_art, Leg_art, # Arterial Tree
  UpBd_vein, SVC, Renal_vein, Splanchnic_vein, Leg_vein, Abd_veins, Thor_IVC, # Venous Tree
  UpBd_cap, Renal_cap, Splanchnic_cap, Leg_cap, # Microcirculation
  Interstitial, # Interstitial Compartment
  Intrathoracic, Abdominal, External, ExternalLBNP, # External Pressures
  ABRafferent, abr_αr, abr_αv, abr_β, abr_para, # Arterial Baroreflex
  CPRafferent, cpr_αr, cpr_αv, # Cardiopulmonary Reflex
  alpha_driver, gravity_driver, lbnp_driver, # Design of Experiments Drivers
  Lung, RespMuscles # Lung Model Breathing ChestWall
  ])

circ_sys = structural_simplify(circ_model)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_model))

# include("InitialUpdate.jl")

"""
Initial Conditions
This section of the code sets the initial conditions for the model. The initial pressures are taken from the vector output of the "initial.jl" file, while the initial interstitial volume is set to a steady-state value. The initial conditions for the reflex arcs are also set here to zero, along with the initial conditions for the sino-atrial node. There will be the same number of initial conditions as there are unknowns, however some of the reflex transfer function ICs are vectorized so the number may not match.
"""

u0 = [
  #### Hemodynamic Pressures
  Asc_A.C.p => x[1],
  BC_A.C.p => x[2],
  UpBd_art.C.p => x[3],
  UpBd_vein.C.p => x[4],
  SVC.C.p => x[5],
  Thor_A.C.p => x[6],
  Abd_A.C.p => x[7],
  Renal_art.C.p => x[8],
  Renal_vein.C.p => x[9],
  Splanchnic_art.C.p => x[10],
  Splanchnic_vein.C.p => x[11],
  Leg_art.C.p => x[12],
  Leg_vein.C.p => x[13],
  Abd_veins.C.p => x[14],
  Thor_IVC.C.p => x[15],
  RA.p => x[16],
  RV.p => x[17], # x[18] is the RV end-systolic pressure
  Pulm_art.C.p => x[19],
  Pulm_vein.C.p => x[20],
  LA.p => x[21],
  LV.p => x[22], # x[23] is the LV end-systolic pressure
  Cor_art.C.p => x[24],
  Cor_vein.C.p => x[25],

  #### Inductance Flows (set to zero for end-diastole)
  BC_A.L.q => 0.0,
  UpBd_art.L.q => 0.0,
  Thor_A.L.q => 0.0,
  Abd_A.L.q => 0.0,
  Renal_art.L.q => 0.0,
  Splanchnic_art.L.q => 0.0,
  Leg_art.L.q => 0.0,
  Cor_art.L.q => 0.0,

  #### Interstitial Compartment
  Interstitial.Qint => 0.0,
  Interstitial.Vint => VintIC,

  #### Reflex Afferents
  CPRafferent.x => x0,
  ABRafferent.x => x0,

  #### Reflex Transfer Functions
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
  cpr_αv.tfend.double_integrator.y => 0.0,

  #### Sino-Atrial Node
  SA.RR_held => RRₙₒₘ,
  SA.ϕ => 0.0,
  SA.Eabr_rv_held => 0.0,
  SA.Eabr_lv_held => 0.0,

  #### Valves
  R_tricuspid.ζ => 1,
  R_pulmonary.ζ => 0,
  R_pulmonary.q => 0,
  R_mitral.ζ => 1,
  R_aortic.ζ => 0,
  R_aortic.q => 0,

  #### Lung Mechanics
  Lung.p_l => p₀,
  Lung.p_tr => p₀,
  Lung.p_b => p₀,
  Lung.p_A => p₀,
  RespMuscles.ϕ => 0.0,
  Lung.pₚₗ => pₚₗ₀
]

"""
SOLVE
This section of the code solves the ODE system. The time span and save interval are set to the values defined in the DOE file.
"""

prob = ODEProblem(circ_sys, u0, tspan)

@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)
# For a discussion on the solver choice, see: https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
# Previously used KenCarp4() which incorporates a stiffness detection and auto-switching algorithm. Tsit5() is the default choice and faster.

"""
Post Solve Calculations
This section calculates the beat-to-beat values, along with some global values (e.g., total blood volume to check for conservation, total peripheral resistance, etc.).
"""

include("PostSolve.jl")

"""
Plots
This section of the code generates some useful plots. Note: all variables in the model can be interrogated.
"""

display(plot(Sol, idxs=[Vtotal],
        label = "TBV",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Total Blood Volume")) # Debugging plot to quickly check volume conservation

#### Direct from Solution Plots

p0a = plot(Sol, idxs=[alpha_driver.α * 180 / π],
             label = "Tilt Angle",
             xlabel = "Time (s)",
             ylabel = "Angle (Degrees)",
             title = "Tilt Angle",
             ylims = (0, 90))
p0b = plot(Sol, idxs=[gravity_driver.g],
             label = "Gravity",
             xlabel = "Time (s)",
             ylabel = "g (m/s^2)",
             title = "Altered-Gravity",
             ylims = (0, 20))
p0c = plot(Sol, idxs=[lbnp_driver.p_lbnp],
             label = "p_LBNP",
             xlabel = "Time (s)",
             ylabel = "LBNP (mmHg)",
             title = "Lower Body Negative Pressure",
             ylims = (-50, 0))

display(plot(p0a, p0b, p0c, layout=(2,2), size=(900,600), suptitle="Design of Experiments"))

p1a = plot(Sol, idxs=[TPR],
             label = "TPR",
             xlabel = "Time (s)",
             ylabel = "TPR (mmHg·s/ml)",
             title = "Microvascular Resistance")
p1b = plot(Sol, idxs=[HR],
             label = "HR",
             xlabel = "Time (s)",
             ylabel = "HR (bpm)",
             title = "Heart Rate")
p1c = plot(Sol, idxs=[Vzpf],
             label = "Vzpf",
             xlabel = "Time (s)",
             ylabel = "Zero-Pressure Volume (ml)",
             title = "Venous Tone")
p1d = plot(Sol, idxs=[RV.Eₘₐₓeff, LV.Eₘₐₓeff],
             label = ["RV" "LV"],
             xlabel = "Time (s)",
             ylabel = "End-Systolic E (mmHg/ml)",
             title = "Ventricular Contractility")

display(plot(p1a, p1b, p1c, p1d, layout=(2,2), size=(900,600), suptitle="Reflex Action"))

p2a = plot(Sol, idxs=[Vtotal, Vvessel, Vinterstitial],
             label = ["Vtotal" "Vvessel" "Vinterstitial"],
             xlabel = "Time (s)",
             ylabel = "Volume (ml)",
             title = "Interstitial Volume")
p2b = plot(Sol, idxs=[LV.pₜₘ, LA.pₜₘ, Asc_A.C.pₜₘ],
             label = ["LV" "LA" "Asc_A"],
             xlabel = "Time (s)",
             ylabel = "Pressure (mmHg)",
             xlims = (50,60),
             title = "Left Heart")
p2c = plot(Sol, idxs=[RA.E, RV.E, LA.E, LV.E],
             label = ["RA" "RV" "LA" "LV"],
             xlabel = "Time (s)",
             ylabel = "E (mmHg/ml)",
             xlims = (50,60),
             title = "Cardiac Elastance")
p2d = plot(Sol, idxs=[Asc_A.in.q],
             label = "Qlvo",
             xlabel = "Time (s)",
             ylabel = "Q (ml/s)",
             xlims = (50,60),
             title = "Left Ventricular Outflow")

display(plot(p2a, p2b, p2c, p2d, layout=(2,2), size=(900,600), suptitle="Hemodynamics"))

#### Beat-to-Beat Plots

p3a = plot(beat_times, [SBP, MAP_integrated, DBP],
             label = ["SBP" "MAP (Integrated)" "DBP"],
             xlabel = "Time (s)",
             ylabel = "Pressure (mmHg)",
             title = "Arterial Pressure",
             ylims = (0, 160))
p3b = plot(beat_times, [CVP],
             label = "CVP",
             xlabel = "Time (s)",
             ylabel = "Pressure (mmHg)",
             title = "Central Venous Pressure")
p3c = plot(beat_times, [CO_Calculated],
             label = "CO",
             xlabel = "Time (s)",
             ylabel = "CO (l/min)",
             title = "Cardiac Output")
p3d = plot(beat_times, [SV],
             label = "SV",
             xlabel = "Time (s)",
             ylabel = "SV (ml)",
             title = "Stroke Volume")

display(plot(p3a, p3b, p3c, p3d, layout=(2,2), size=(900,600), suptitle="Beat-to-Beat Trends"))

display(plot(beat_times, [(UpBd_art_Vmean + UpBd_vein_Vmean),
        (Renal_art_Vmean + Renal_vein_Vmean),
        (Splanchnic_art_Vmean + Splanchnic_vein_Vmean),
        (Leg_art_Vmean + Leg_vein_Vmean),
        (Cor_art_Vmean + Cor_vein_Vmean)],
        label = ["Upper Body" "Renal" "Splanchnic" "Leg" "Coronary"],
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Average Branch Volumes"))

#### Pulmonary and Respiratory Plots

p4a = plot(Sol, idxs=[RespMuscles.out.p], xlims = (0, 250),
        label = "pₘᵤₛ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Respiratory Muscle Pressure")

p4b = plot(Sol, idxs=[Lung.pₚₗ], xlims = (0, 250),
        label = "pₚₗ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Pleural Pressure")

p4c = plot(Sol, idxs=[Lung.p_A], xlims = (0, 250),
        label = "p_A",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Alveolar Pressure")

p4d = plot(Sol, idxs=[Lung.Vrᵢₙ], xlims = (0, 250),
        label = "Air Flow",
        xlabel = "Time (s)",
        ylabel = "Flow (ml/s)",
        title = "Air Flow")

p4e = plot(Sol, idxs=[Lung.V_A + Lung.V_D], xlims = (0, 250),
        label = "V_L",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Lung Volume")

p4f = plot(Sol, idxs=[Lung.V_A], xlims = (0, 250),
        label = "V_A",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Alveolar Volume")

p4g = plot(Sol, idxs=[Lung.V_D], xlims = (0, 250),
        label = "V_D",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Dead Space Volume")

p4h = plot(Sol, idxs=[Intrathoracic.pth.p], xlims = (0, 250),
      label = "pₜₕ",
      xlabel = "Time (s)",
      ylabel = "Pressure (mmHg)",
      title = "Intrathoracic Pressure")

display(plot(p4a,p4b,p4c,p4d,p4e,p4f,p4g,p4h, layout=(4,2), size=(900,600), suptitle="Lungs"))

"""
Save Outputs
Uncomment the following lines to save the outputs to a CSV file.
"""

# using CSV
# using DataFrames

# df_cycle = DataFrame(
#     beat_times = beat_times,
#     SBP = SBP,
#     DBP = DBP,
#     MAP = MAP_integrated,
#     PP = PP,
#     CVP = CVP,
#     SV = SV,
#     CO = CO_Calculated,
#     HR = HR_beats,
#     RR_int = RR_beats,
#     TPR = TPR_beats,
#     EF = EF
# ) # Create a DataFrame with beat-by-beat cardiovascular metrics

# CSV.write("output.csv", df_cycle, header=true, delim=',', append=false) # Write beat-by-beat metrics to a CSV file


# all_para_names = names(ModelParams, all=true) # Extract all parameter names from the ModelParameters.jl file
# params_syms = filter(n -> isdefined(ModelParams, n) && getfield(ModelParams, n) isa Number, all_para_names) # Filter for numerical (scalar) parameters
# parameter_names  = String.(params_syms) # Convert symbols to strings for DataFrame labeling
# parameter_values = [ getfield(ModelParams, s) for s in params_syms ] # Retrieve corresponding parameter values
# df_parameters = DataFrame(Parameter = parameter_names, Value = parameter_values) # Construct a DataFrame of parameter names and values

# CSV.write("model_parameters.csv", df_parameters) # Write model parameters to CSV

end