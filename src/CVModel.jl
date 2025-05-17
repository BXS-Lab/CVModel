"""
Cardiovascular Model
Version 3.3.0 (May 14th, 2025)
BXS Lab, UC Davis; Authors: RS Whittle, AJ Kondoor, HS Vellore
Contact info: Dr. Rich Whittle – Department of Mechanical and Aerospace Engineering, UC Davis, One Shields Ave, Davis CA 95616 (rswhittle@ucdavis.edu)
This model is a simulation of the human cardiovascular system, including a four chamber heart, arteries, veins, and microcirculation. It includes reflex control for the arterial baroreflex (ABR) and cardiopulmonary reflex (CPR), as well as hydrostatic effects, interstitial fluid flow, and external tissue pressures. The model is designed to simulate various physiological scenarios, including a tilt angle protocol, altered-gravity environment, and lower body negative pressure (LBNP) protocol. The underlying equations are based on the work of Heldt (2004), Zamanian (2007), Diaz Artiles (2015), Albanese (2016), and Whittle (2023). The model is implemented using the ModelingToolkit.jl package in Julia.
"""
### TODO: CV Model:            (1) Vertebral Plexus (2) Dynamic ICP
### TODO: Pulmonary Mechanics: (1) Bring Intrathoracic Pressure into Lung Model
### TODO: Simulation:          (1) Exercise Model, (2) Other blood parameters (e.g., pH etc.), (https://github.com/baillielab/oxygen_delivery/blob/master/oxygen_delivery.py) (3) Altitude, pressure, temperature driver; water vapor etc. (4) Update ICs (negative compliance alters volume)
### TODO: Code:                (1) Fix the whole model params thing
module CVModel
display("Cardiovascular Model v3.3.0 (May 14th, 2025) - BXS Lab, UC Davis")

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
using NLsolve

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
@named Pulm_art = Artery(C=Cpa, R=Rpa, V₀=v0pa, has_hydrostatic=false, has_tissue=false, has_inertia=false, is_pulmonary=true)
@named Pulm_cap = StarlingResistor(R=Rpc, h=h_Lungs)
@named Pulm_vein = Vein(C=Cpv, R=Rpv, V₀=v0pv, has_hydrostatic=false, has_tissue=false)

#### Left Heart
@named LA = HeldtChamber(V₀=v0_la, Eₘᵢₙ=Ed_la, Eₘₐₓ=Ees_la, τₑₛ=τₐₛ, inP=true)
@named R_mitral = MynardValve_Atrioventricular(Ann=Ann_mv, Kvc=Kvc_mv, Kvo=Kvo_mv)
@named LV = HeldtChamber(V₀=v0_lv, Eₘᵢₙ=Ed_lv, Eₘₐₓ=Ees_lv, Elimit = Elimit_lv, τₑₛ=τᵥₛ, inP=true, has_abr=true)
@named R_aortic = MynardValve_SemiLunar(Leff=Leff_av, Ann=Ann_av, Kvc=Kvc_av, Kvo=Kvo_av)

#### Coronary Circulation
@named Cor_art = Artery(C=Cca, R=Rca, V₀=v0ca, has_hydrostatic=false, has_tissue=false, is_heart=true, Vₜ=Vₜ_heart, RQ=RQ₀)
@named Cor_cap = VarResistor()
@named Cor_vein = Vein(C=Ccv, R=Rcv, V₀=v0cv, has_hydrostatic=false, has_tissue=false)

#### Dynamic Heart Power
@named HeartP = HeartPower()

#### Arterial Circulation
@named Asc_A = Artery(C=C_Asc_A, R=R_Asc_A, V₀=v0_Asc_A, h=h_Asc_A, has_tissue=false, has_inertia=false)
@named BC_A = Artery(C=C_BC_A, R=R_BC_A, V₀=v0_BC_A, h=h_BC_A, has_tissue=false, L=L_BC_A)
@named UpBd_art = Artery(C=C_UpBd_art, R=R_UpBd_art, V₀=v0_UpBd_art, h=h_UpBd_art, con=con_UpBd_art, rad=rad_UB, L=L_UpBd_art, has_gasexchange=true, Vₜ=Vₜ_ub, MO₂=MO₂_ub, RQ=RQ₀)
@named Thor_A = Artery(C=C_Thor_A, R=R_Thor_A, V₀=v0_Thor_A, h=h_Thor_A, has_tissue=false, L=L_Thor_A)
@named Abd_A = Artery(C=C_Abd_A, R=R_Abd_A, V₀=v0_Abd_A, h=h_Abd_A, rad=rad_Abd, L=L_Abd_A)
@named Renal_art = Artery(C=C_Renal_art, R=R_Renal_art, V₀=v0_Renal_art, h=h_Renal_art, rad=rad_Abd, L=L_Renal_art, has_gasexchange=true, Vₜ=Vₜ_renal, MO₂=MO₂_renal, RQ=RQ₀)
@named Splanchnic_art = Artery(C=C_Splanchnic_art, R=R_Splanchnic_art, V₀=v0_Splanchnic_art, h=h_Splanchnic_art, rad=rad_Abd, L=L_Splanchnic_art, has_gasexchange=true, Vₜ=Vₜ_splanchnic, MO₂=MO₂_splanchnic, RQ=RQ₀)
@named Leg_art = Artery(C=C_Leg_art, R=R_Leg_art, V₀=v0_Leg_art, h=h_Leg_art, con=con_Leg_art, rad=rad_Leg, L=L_Leg_art, has_gasexchange=true, Vₜ=Vₜ_legs, MO₂=MO₂_legs, RQ=RQ₀)
@named CommonCarotid = Artery(C=C_CCA, R=R_CCA, V₀=v0_CCA, h=h_CCA, rad=rad_Neck, L=L_CCA)
@named Head_art = Artery(C=C_Head_art, R=R_Head_art, V₀=v0_Head_art, h=h_Head_art, rad=rad_Head, L=L_Head_art, has_gasexchange=true, Vₜ=Vₜ_brain, MO₂=MO₂_brain, RQ=RQ₀)

#### Venous Circulation
@named UpBd_vein = Vein(C=C_UpBd_vein, R=R_UpBd_vein, V₀=v0_UpBd_vein, h=h_UpBd_vein, con=con_UpBd_vein, has_valve=true, has_reflex=true, rad=rad_UB)
@named SVC = Vein(C=C_SVC, R=R_SVC, V₀=v0_SVC, h=h_SVC, has_tissue=false)
@named Renal_vein = Vein(C=C_Renal_vein, R=R_Renal_vein, V₀=v0_Renal_vein, h=h_Renal_vein, has_reflex=true, rad=rad_Abd, V_min=vₘᵢₙ_Renal_vein)
@named Splanchnic_vein = Vein(C=C_Splanchnic_vein, R=R_Splanchnic_vein, V₀=v0_Splanchnic_vein, h=h_Splanchnic_vein, is_nonlinear=true, Flow_div = Flow_Splanchnic_vein, V_max=vM_Splanchnic_vein, has_reflex=true, rad=rad_Abd)
@named Leg_vein = Vein(C=C_Leg_vein, R=R_Leg_vein, V₀=v0_Leg_vein, h=h_Leg_vein, con=con_Leg_vein, has_valve=true, is_nonlinear=true, Flow_div = Flow_Leg_vein, V_max=vM_Leg_vein, has_reflex=true, rad=rad_Leg)
@named Abd_veins = Vein(C=C_Abd_veins, R=R_Abd_veins, V₀=v0_Abd_veins, h=h_Abd_veins, is_nonlinear=true, Flow_div = Flow_Abd_veins, V_max=vM_Abd_vein, rad=rad_Abd)
@named Thor_IVC = Vein(C=C_Thor_IVC, R=R_Thor_IVC, V₀=v0_Thor_IVC, h=h_Thor_IVC, has_tissue=false)
@named Head_veins = Vein(C=C_Head_veins, R=R_Head_veins, V₀=v0_Head_veins, h=h_Head_veins, rad=rad_Head)
@named Jugular_vein = Vein(C=C_Jugular_vein, R=R_Jugular_vein, V₀=v0_Jugular_vein, h=h_Jugular_vein, rad=rad_Neck, has_valve=true)

#### Vascular Junctions (necessary for blood gas)
@named Asc_A_Junc = Junction3()
@named BC_A_Junc = Junction2()
@named Abd_A_Junc = Junction3()
@named Head_veins_Junc = Junction2()

# @named VP = VertebralPlexus(R=Rᵥₚ, h=h_Jugular_vein)

#### Microcirculation
@named Head_cap = VarConductor()
@named UpBd_cap = VarResistor()
@named Renal_cap = VarResistor()
@named Splanchnic_cap = VarResistor()
@named Leg_cap = VarResistor()

#### Interstitial Compartment
@named Interstitial = InterstitialCompartment(Vmtilt=Vmax_tilt, Vmlbnp=Vmax_lbnp, τ=τint)

#### External Pressures
@named External = ExternalPressureUB(p_ext=p₀)
@named Intracranial = IntracranialPressure(p_icp=p_icp)
@named Intrathoracic = IntrathoracicPressure()
@named Abdominal = IntraAbdominalPressure(p_abd=p_abd)
@named ExternalLBNP = ExternalPressureLB(p_ext=p₀)

#### Lung model
@named RespMuscles = RespiratoryMuscles()
@named Lungs = Lung()
@named LungGE = LungGasExchange()
@named TV = TidalVolume()

################
# Reflex Arcs
################

#### Autoregulation
@named BrainAutoreg = CerebralAutoregulation()
@named HeartAutoreg = Autoregulation(_gjO₂=ghO₂, _CvjO₂n=CvhO₂n, _τO₂=τO₂, _PaCO₂n=PaCO₂n, _kjCO₂=khCO₂, _τCO₂=τCO₂)
@named UBMuscleAutoreg = Autoregulation(_gjO₂=gmO₂, _CvjO₂n=CvmO₂n, _τO₂=τO₂, _PaCO₂n=PaCO₂n, _kjCO₂=kmCO₂, _τCO₂=τCO₂)
@named LBMuscleAutoreg = Autoregulation(_gjO₂=gmO₂, _CvjO₂n=CvmO₂n, _τO₂=τO₂, _PaCO₂n=PaCO₂n, _kjCO₂=kmCO₂, _τCO₂=τCO₂)

#### Arterial Baroreflex
@named ABR = AfferentBaroreflex()

#### Cardiopulmonary Reflex
@named CPR = AfferentCPR()

#### Peripheral Chemoreceptors
@named PeripheralChemo = PeripheralChemoreceptors()

#### Lung Stretch Receptors
@named LungStretchReceptors = LungStretch()

#### CNS Ischemic Response
@named IschArterioles = IschemicResponse(_χₛⱼ=χₛₚ, _PaO₂ₛⱼn=PO₂nₛₚ, _kiscₛⱼ=kiscₛₚ, _τisc=τisc, _PaCO₂n=PaCO₂n, _gccₛⱼ=gccₛₚ, _τcc=τcc, _θₛⱼn=θₛₚₙ)
@named IschVeins = IschemicResponse(_χₛⱼ=χₛᵥ, _PaO₂ₛⱼn=PO₂nₛᵥ, _kiscₛⱼ=kiscₛᵥ, _τisc=τisc, _PaCO₂n=PaCO₂n, _gccₛⱼ=gccₛᵥ, _τcc=τcc, _θₛⱼn=θₛᵥₙ)
@named IschHeart = IschemicResponse(_χₛⱼ=χₛₕ, _PaO₂ₛⱼn=PO₂nₛₕ, _kiscₛⱼ=kiscₛₕ, _τisc=τisc, _PaCO₂n=PaCO₂n, _gccₛⱼ=gccₛₕ, _τcc=τcc, _θₛⱼn=θₛₕₙ)

#### Efferent Pathways
@named Efferent = EfferentPathways()

#### Effectors
@named ERH = Effectors(Gain=GEmaxrv, delay=DEmaxrv, time=τEmaxrv, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named ELH = Effectors(Gain=GEmaxlv, delay=DEmaxlv, time=τEmaxlv, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named ERR = EffectorsRR(Gainₛ=GTs, Gainᵥ=GTv, delayₛ=DTs, delayᵥ=DTv, timeₛ=τTs, timeᵥ=τTv, min=fesₘᵢₙ, delay_order=reflex_delay_order)

@named EResistance_UpBd = Effectors(Gain=GResistance_UpBd, delay=DResistance, time=τResistance, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EResistance_Renal = Effectors(Gain=GResistance_Renal, delay=DResistance, time=τResistance, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EResistance_Splanchnic = Effectors(Gain=GResistance_Splanchnic, delay=DResistance, time=τResistance, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EResistance_Leg = Effectors(Gain=GResistance_Leg, delay=DResistance, time=τResistance, min=fesₘᵢₙ, delay_order=reflex_delay_order)

@named EVtone_UpBd = Effectors(Gain=GVtone_UpBd, delay=DVtone, time=τVtone, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EVtone_Renal = Effectors(Gain=GVtone_Renal, delay=DVtone, time=τVtone, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EVtone_Splanchnic = Effectors(Gain=GVtone_Splanchnic, delay=DVtone, time=τVtone, min=fesₘᵢₙ, delay_order=reflex_delay_order)
@named EVtone_Leg = Effectors(Gain=GVtone_Leg, delay=DVtone, time=τVtone, min=fesₘᵢₙ, delay_order=reflex_delay_order)

#### Pulmonary Reflexes
@named CentralResp = Chemoreceptors(Delay=Dc, Gain_A=G_cA, Gain_f=G_cf, set_point=PaCO₂n, time_A=τ_cA, time_f=τ_cf, delay_order=reflex_delay_order)
@named PeripheralResp = Chemoreceptors(Delay=Dp, Gain_A=G_pA, Gain_f=G_pf, set_point=fapc_set, time_A=τ_pA, time_f=τ_pf, delay_order=reflex_delay_order)

#### Design of Experiments
# These are the drivers to implement the protocols defined in the "doe.jl" file.
@named alpha_driver = Alpha()
@named gravity_driver = Gravity()
@named lbnp_driver = LBNP()

"""
Structural Connections
This section of code connects the instanced compartments together to form the cardiovascular system. The hemodynamic connections are made using the "connect" function, which connects pins linked in pressure and flow via Kirchhoff's laws. The DOE integrations and reflex connections are just included as direct signal connections to the respective compartments.
"""

circ_eqs = [
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
  LV.τ ~ SA.RR_held,

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
  connect(Asc_A.out, Asc_A_Junc.in),
  connect(Asc_A_Junc.out1, BC_A.in),
  connect(Asc_A_Junc.out2, Thor_A.in),
  connect(Asc_A_Junc.out3, Cor_art.in),
  connect(BC_A.out, BC_A_Junc.in),
  connect(BC_A_Junc.out1, UpBd_art.in),
  connect(BC_A_Junc.out2, CommonCarotid.in),
  connect(Thor_A.out, Abd_A.in),
  connect(Abd_A.out, Abd_A_Junc.in),
  connect(Abd_A_Junc.out1, Renal_art.in),
  connect(Abd_A_Junc.out2, Splanchnic_art.in),
  connect(Abd_A_Junc.out3, Leg_art.in),

  #### Head and Neck Circulation
  connect(CommonCarotid.out, Head_art.in),
  connect(Head_art.out, Head_cap.in),
  connect(Head_cap.out, Head_veins.in),
#   connect(Head_veins.out, Head_veins_Junc.in),
#   connect(Head_veins_Junc.out1, Jugular_vein.in),
#   connect(Head_veins_Junc.out2, VP.in),
  connect(Head_veins.out, Jugular_vein.in),

  #### Coronary Circulation
  connect(Cor_art.out, Cor_cap.in),
  connect(Cor_cap.out, Cor_vein.in),

  #### Heart Power
  HeartP.Plv ~ LV.pₜₘ,
  HeartP.Prv ~ RV.pₜₘ,
  HeartP.DVlv ~ (LV.in.q + LV.out.q),
  HeartP.DVrv ~ (RV.in.q + RV.out.q),
  HeartP.MO₂dyn ~ Cor_art.C.MO₂dyn,

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
#   connect(Jugular_vein.out, VP.out, UpBd_vein.out, SVC.in),
  connect(Jugular_vein.out, UpBd_vein.out, SVC.in),
  connect(Abd_veins.in, Renal_vein.out, Splanchnic_vein.out, Leg_vein.out),
  connect(Abd_veins.out, Thor_IVC.in),
  connect(Thor_IVC.out, SVC.out, Cor_vein.out, RA.in),

  #### External Pressures
  connect(Intrathoracic.pth, Asc_A.ep, BC_A.ep, Thor_A.ep, SVC.ep, Thor_IVC.ep, RA.ep, RV.ep, Pulm_art.ep, Pulm_vein.ep, LA.ep, LV.ep, Cor_art.ep, Cor_vein.ep),
  connect(Abdominal.pabd, Abd_A.ep, Renal_art.ep, Splanchnic_art.ep, Renal_vein.ep, Splanchnic_vein.ep, Abd_veins.ep),
  connect(External.pext, UpBd_art.ep, UpBd_vein.ep, CommonCarotid.ep, Jugular_vein.ep, Lungs.in, RespMuscles.in),
  connect(ExternalLBNP.pext, Leg_art.ep, Leg_vein.ep),
  connect(Intracranial.picp, Head_art.ep, Head_veins.ep),

  #### Interstitial Connections (Direct Connections)
  Splanchnic_vein.C.qint ~ Interstitial.Qint,
  Leg_vein.C.qint ~ Interstitial.Qint,
  Abd_veins.C.qint ~ Interstitial.Qint,

  #### Pulmonary Connections
  connect(Lungs.chestwall, RespMuscles.out),
  TV.breath_trigger ~ RespMuscles.breath_trigger,
  TV.breath_trigger2 ~ RespMuscles.breath_trigger2,
  TV.V_L ~ Lungs.V_L,

  #### Pulmonary to Cardiovascular Pressure Connections
  Pulm_cap.pₐₗᵥ ~ Lungs.p_A,
  Intrathoracic.pₚₗ ~ Lungs.pₚₗ,

  #### Pulmary Gas Exchenge
  LungGE.Vrᵢₙ ~ Lungs.Vrᵢₙ,
  LungGE.Vr_A ~ Lungs.Vr_A,
  LungGE.V_D ~ Lungs.V_D,
  LungGE.V_A ~ Lungs.V_A,
  LungGE.qpa ~ Pulm_art.q,
  LungGE.Vpa ~ Pulm_art.C.V,
  LungGE.cvO₂ ~ Pulm_art.C.cO₂,
  LungGE.cvCO₂ ~ Pulm_art.C.cCO₂,
  LungGE.caO₂ ~ Pulm_art.C.caO₂,
  LungGE.caCO₂ ~ Pulm_art.C.caCO₂,

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
  CommonCarotid.α ~ alpha_driver.α,
  Head_art.α ~ alpha_driver.α,
  Head_veins.α ~ alpha_driver.α,
  Jugular_vein.α ~ alpha_driver.α,
#   VP.α ~ alpha_driver.α,
  Interstitial.α ~ alpha_driver.α,
  Intrathoracic.α ~ alpha_driver.α,
  Pulm_art.α ~ alpha_driver.α,
  Pulm_vein.α ~ alpha_driver.α,
  Pulm_cap.α ~ alpha_driver.α,
  Lungs.α ~ alpha_driver.α,

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
  CommonCarotid.g ~ gravity_driver.g,
  Head_art.g ~ gravity_driver.g,
  Head_veins.g ~ gravity_driver.g,
  Jugular_vein.g ~ gravity_driver.g,
#   VP.g ~ gravity_driver.g,
  Interstitial.g ~ gravity_driver.g,
  Intrathoracic.g ~ gravity_driver.g,
  Pulm_art.g ~ gravity_driver.g,
  Pulm_vein.g ~ gravity_driver.g,
  Pulm_cap.g ~ gravity_driver.g,
  Lungs.g ~ gravity_driver.g,

  #### LBNP Equations (Direct Connections)
  ExternalLBNP.p_lbnp ~ lbnp_driver.p_lbnp,
  Interstitial.p_lbnp ~ lbnp_driver.p_lbnp,

  ################
  # Reflex Arcs
  ################

  #### Autoregulation
  BrainAutoreg.uCvbO₂ ~ Head_veins.cO₂,
  HeartAutoreg.uCvjO₂ ~ Cor_vein.cO₂,
  UBMuscleAutoreg.uCvjO₂ ~ UpBd_vein.cO₂,
  LBMuscleAutoreg.uCvjO₂ ~ Leg_vein.cO₂,
  BrainAutoreg.uPaCO₂ ~ LungGE.paCO₂,
  HeartAutoreg.uPaCO₂ ~ LungGE.paCO₂,
  UBMuscleAutoreg.uPaCO₂ ~ LungGE.paCO₂,
  LBMuscleAutoreg.uPaCO₂ ~ LungGE.paCO₂,

  #### Afferent: Arterial Baroreflex
  ABR.pb ~ (Asc_A.C.pₜₘ+(BC_A.C.pₜₘ + (ρ_b*gravity_driver.g*(h_cs/100)*sin(alpha_driver.α)*Pa2mmHg)))/2.0,

  #### Afferent: Cardiopulmonary Reflex
  CPR.pr ~ RA.pₜₘ,

  #### Afferent: Peripheral Chemoreceptors
  PeripheralChemo.uSaO₂ ~ LungGE.SaO₂,
  PeripheralChemo.ucaCO₂ ~ LungGE.caCO₂,

  #### Afferent: Lung Stretch Receptors
  TV.VT ~ LungStretchReceptors.VT,

  #### Afferent: CNS Ischemic Response
  IschArterioles.uPaO₂ ~ LungGE.paO₂,
  IschArterioles.uPaCO₂ ~ LungGE.paCO₂,
  IschVeins.uPaO₂ ~ LungGE.paO₂,
  IschVeins.uPaCO₂ ~ LungGE.paCO₂,
  IschHeart.uPaO₂ ~ LungGE.paO₂,
  IschHeart.uPaCO₂ ~ LungGE.paCO₂,

  #### Efferent Pathways
  Efferent.fab ~ ABR.fab,
  Efferent.fcpr ~ CPR.fcpr,
  Efferent.fapc ~ PeripheralChemo.fapc,
  Efferent.fasr ~ LungStretchReceptors.fasr,
  Efferent.θₛₕ ~ IschHeart.θₛⱼ,
  Efferent.θₛₚ ~ IschArterioles.θₛⱼ,
  Efferent.θₛᵥ ~ IschVeins.θₛⱼ,

  #### Effectors
  ERH.u ~ Efferent.fₛₕ,
  ELH.u ~ Efferent.fₛₕ,
  ERR.uₛ ~ Efferent.fₛₕ,
  ERR.uᵥ ~ Efferent.fᵥ,

  EResistance_UpBd.u ~ Efferent.fₛₚ,
  EResistance_Renal.u ~ Efferent.fₛₚ,
  EResistance_Splanchnic.u ~ Efferent.fₛₚ,
  EResistance_Leg.u ~ Efferent.fₛₚ,

  EVtone_UpBd.u ~ Efferent.fₛᵥ,
  EVtone_Renal.u ~ Efferent.fₛᵥ,
  EVtone_Splanchnic.u ~ Efferent.fₛᵥ,
  EVtone_Leg.u ~ Efferent.fₛᵥ,

  #### Effectors: Arteriole Resistance
  Cor_cap.R ~ Rcc * (1 + HeartAutoreg.xjCO₂) /(1 + HeartAutoreg.xjO₂),
  Head_cap.G ~ Gbpn * (1 + BrainAutoreg.xbO₂ + BrainAutoreg.xbCO₂),

  UpBd_cap.R ~ (R_UpBd_cap + EResistance_UpBd.Δσ) * (1 + UBMuscleAutoreg.xjCO₂) /(1 + UBMuscleAutoreg.xjO₂),
  Renal_cap.R ~ R_Renal_cap + EResistance_Renal.Δσ,
  Splanchnic_cap.R ~ R_Splanchnic_cap + EResistance_Splanchnic.Δσ,
  Leg_cap.R ~ (R_Leg_cap + EResistance_Leg.Δσ) * (1 + LBMuscleAutoreg.xjCO₂) /(1 + LBMuscleAutoreg.xjO₂),

  #### Effectors: Venous Tone
  UpBd_vein.ΔV ~ EVtone_UpBd.Δσ,
  Renal_vein.ΔV ~ EVtone_Renal.Δσ,
  Splanchnic_vein.ΔV ~ EVtone_Splanchnic.Δσ,
  Leg_vein.ΔV ~ EVtone_Leg.Δσ,

  #### Effectors: Ventricular Contractility
  SA.Eabr_rv ~ ERH.Δσ,
  SA.Eabr_lv ~ ELH.Δσ,

  #### Effectors: Heart Rate
  SA.RRabr ~ ERR.ΔT,

  #### Pulmonary Reflexes
  CentralResp.u ~ LungGE.paCO₂,
  PeripheralResp.u ~ PeripheralChemo.fapc,

  RespMuscles.RespRate_chemo ~ (CentralResp.y_f + PeripheralResp.y_f),
  RespMuscles.p_chemo ~ (CentralResp.y_A + PeripheralResp.y_A),
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
  CommonCarotid, Head_art, Head_veins, Jugular_vein, # Head and Neck Circulation
#   VP, # Vertebral Plexus
  Asc_A_Junc, BC_A_Junc, Abd_A_Junc, # Arterial Junctions
#   Head_veins_Junc, # Venous Junctions
  UpBd_cap, Renal_cap, Splanchnic_cap, Leg_cap, Head_cap, # Microcirculation
  Interstitial, # Interstitial Compartment
  Intrathoracic, Abdominal, External, ExternalLBNP, Intracranial, # External Pressures
  alpha_driver, gravity_driver, lbnp_driver, # Design of Experiments Drivers
  Lungs, RespMuscles, LungGE, TV, # Lung Model Breathing ChestWall
  CentralResp, PeripheralResp,
  PeripheralChemo, LungStretchReceptors, # Respiratory Control
  BrainAutoreg, HeartAutoreg, UBMuscleAutoreg, LBMuscleAutoreg,# Autoregulation
  HeartP, # Heart Power
  ABR, CPR, # Afferent Baroreflex
  IschArterioles, IschVeins, IschHeart, # CNS Ischemic Response
  Efferent, # Efferent Pathways
  ERH, ELH, ERR, # Effectors
  EResistance_UpBd, EResistance_Renal, EResistance_Splanchnic, EResistance_Leg, # Effectors: Arteriole Resistance
  EVtone_UpBd, EVtone_Renal, EVtone_Splanchnic, EVtone_Leg, # Effectors: Venous Tone
  ])

circ_sys = structural_simplify(circ_model)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
equations(expand(circ_sys))
unknowns(circ_sys)
equations(expand(circ_model))
# unknown_list = collect(unknowns(circ_sys))
# println("Number of unknowns: ", length(unknown_list))
# for (i, u) in enumerate(unknown_list)
#     println("unknown $(i): ", u)
# end

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
  CommonCarotid.C.p => x[16],
  Head_art.C.p => x[17],
  Head_veins.C.p => x[18],
  Jugular_vein.C.p => x[19],
  RA.p => x[20],
  RV.p => x[21], # x[22] is the RV end-systolic pressure
  Pulm_art.C.p => x[23],
  Pulm_vein.C.p => x[24],
  LA.p => x[25],
  LV.p => x[26], # x[27] is the LV end-systolic pressure
  Cor_art.C.p => x[28],
  Cor_vein.C.p => x[29],

  #### Inductance Flows (set to zero for end-diastole)
  BC_A.L.q => 0.0,
  UpBd_art.L.q => 0.0,
  Thor_A.L.q => 0.0,
  Abd_A.L.q => 0.0,
  Renal_art.L.q => 0.0,
  Splanchnic_art.L.q => 0.0,
  Leg_art.L.q => 0.0,
  CommonCarotid.L.q => 0.0,
  Head_art.L.q => 0.0,
  Cor_art.L.q => 0.0,

  #### Interstitial Compartment
  Interstitial.Qint => 0.0,
  Interstitial.Vint => VintIC,

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
  RespMuscles.ϕ => 0.0,
  Lungs.p_l => p₀,
  Lungs.p_tr => p₀,
  Lungs.p_b => p₀,
  Lungs.p_A => p₀,
  Lungs.pₚₗ => pₚₗ₀,
  TV.VL_min => 2500,
  TV.VL_max => 3000,

  #### Blood Gas O₂ Arterial
#   Pulm_art.C.out.cO₂ => 0.0,
  Pulm_vein.C.out.cO₂ => 0.2,
  LA.out.cO₂ => 0.2,
  LV.out.cO₂ => 0.2,
  Cor_art.C.out.cO₂ => 0.152, # After exchange
  Asc_A.C.out.cO₂ => 0.2,
  BC_A.C.out.cO₂ => 0.2,
  UpBd_art.C.out.cO₂ => 0.152, # After exchange
  Thor_A.C.out.cO₂ => 0.2,
  Abd_A.C.out.cO₂ => 0.2,
  Renal_art.C.out.cO₂ => 0.152, # After exchange
  Splanchnic_art.C.out.cO₂ => 0.152, # After exchange
  Leg_art.C.out.cO₂ => 0.152, # After exchange
  CommonCarotid.C.out.cO₂ => 0.2,
  Head_art.C.out.cO₂ => 0.152, # After exchange
  #### Blood Gas O₂ Venous
  Cor_vein.C.out.cO₂ => 0.152,
  UpBd_vein.C.out.cO₂ => 0.152,
  SVC.C.out.cO₂ => 0.152,
  Renal_vein.C.out.cO₂ => 0.152,
  Splanchnic_vein.C.out.cO₂ => 0.152,
  Leg_vein.C.out.cO₂ => 0.152,
  Abd_veins.C.out.cO₂ => 0.152,
  Thor_IVC.C.out.cO₂ => 0.152,
  Head_veins.C.out.cO₂ => 0.152,
  Jugular_vein.C.out.cO₂ => 0.152,
  RA.out.cO₂ => 0.152,
  RV.out.cO₂ => 0.152,

  #### Blood Gas CO₂ Arterial
#   Pulm_art.C.out.cCO₂ => 0.0,
  Pulm_vein.C.out.cCO₂ => 0.485,
  LA.out.cCO₂ => 0.485,
  LV.out.cCO₂ => 0.485,
  Cor_art.C.out.cCO₂ => 0.529, # After exchange
  Asc_A.C.out.cCO₂ => 0.485,
  BC_A.C.out.cCO₂ => 0.485,
  UpBd_art.C.out.cCO₂ => 0.529, # After exchange
  Thor_A.C.out.cCO₂ => 0.485,
  Abd_A.C.out.cCO₂ => 0.485,
  Renal_art.C.out.cCO₂ => 0.529, # After exchange
  Splanchnic_art.C.out.cCO₂ => 0.529, # After exchange
  Leg_art.C.out.cCO₂ => 0.529, # After exchange
  CommonCarotid.C.out.cCO₂ => 0.485,
  Head_art.C.out.cCO₂ => 0.529, # After exchange
  #### Blood Gas CO₂ Venous
  Cor_vein.C.out.cCO₂ => 0.529,
  UpBd_vein.C.out.cCO₂ => 0.529,
  SVC.C.out.cCO₂ => 0.529,
  Renal_vein.C.out.cCO₂ => 0.529,
  Splanchnic_vein.C.out.cCO₂ => 0.529,
  Leg_vein.C.out.cCO₂ => 0.529,
  Abd_veins.C.out.cCO₂ => 0.529,
  Thor_IVC.C.out.cCO₂ => 0.529,
  Head_veins.C.out.cCO₂ => 0.529,
  Jugular_vein.C.out.cCO₂ => 0.529,
  RA.out.cCO₂ => 0.529,
  RV.out.cCO₂ => 0.529,

  #### Lung Fractional Concentrations
  LungGE.FDO₂ => FIO₂,
  LungGE.FDCO₂ => FICO₂,
  LungGE.FAO₂ => FIO₂,
  LungGE.FACO₂ => FICO₂,

  #### Peripheral Chemoreceptors
  PeripheralChemo.ϕCO₂dyn => 0.0,
  PeripheralChemo.ϕc => 0.0,

  #### Lung Stretch Receptors
  LungStretchReceptors.fasr => 0.0,

  #### Respiratory Control
  CentralResp.delay.x => reflex_delay_init,
  CentralResp.y_A => 0.0,
  CentralResp.y_f => 0.0,
  PeripheralResp.delay.x => reflex_delay_init,
  PeripheralResp.y_A => 0.0,
  PeripheralResp.y_f => 0.0,

  #### Central Pattern Generator
  RespMuscles.BreathInt_held => 60 / RespRateₙₒₘ,
  RespMuscles.p_held => p_musmin,

  #### Autoregulation
  BrainAutoreg.xbO₂ => 0.0,
  BrainAutoreg.xbCO₂ => 0.0,
  HeartAutoreg.xjO₂ => 0.0,
  HeartAutoreg.xjCO₂ => 0.0,
  UBMuscleAutoreg.xjO₂ => 0.0,
  UBMuscleAutoreg.xjCO₂ => 0.0,
  LBMuscleAutoreg.xjO₂ => 0.0,
  LBMuscleAutoreg.xjCO₂ => 0.0,

  #### Heart Dynamic Consumption
  HeartP.Wh => Whₙₒₘ,

  #### Afferent Baroreflex & CPR
  ABR.P => Pn,
  CPR.P => Prn,

  #### Ischemic Response
  IschArterioles.ΔΘO₂ₛⱼ => 0.0,
  IschArterioles.ΔΘCO₂ₛⱼ => 0.0,
  IschVeins.ΔΘO₂ₛⱼ => 0.0,
  IschVeins.ΔΘCO₂ₛⱼ => 0.0,
  IschHeart.ΔΘO₂ₛⱼ => 0.0,
  IschHeart.ΔΘCO₂ₛⱼ => 0.0,

  #### Effectors
  ERH.Δσ => 0.0,
  ERH.d.x => reflex_delay_init,
  ELH.Δσ => 0.0,
  ELH.d.x => reflex_delay_init,
  ERR.ΔT => 0.0,
  ERR.dₛ.x => reflex_delay_init,
  ERR.dᵥ.x => reflex_delay_init,

  EResistance_UpBd.Δσ => 0.0,
  EResistance_UpBd.d.x => reflex_delay_init,
  EResistance_Renal.Δσ => 0.0,
  EResistance_Renal.d.x => reflex_delay_init,
  EResistance_Splanchnic.Δσ => 0.0,
  EResistance_Splanchnic.d.x => reflex_delay_init,
  EResistance_Leg.Δσ => 0.0,
  EResistance_Leg.d.x => reflex_delay_init,

  EVtone_UpBd.Δσ => 0.0,
  EVtone_UpBd.d.x => reflex_delay_init,
  EVtone_Renal.Δσ => 0.0,
  EVtone_Renal.d.x => reflex_delay_init,
  EVtone_Splanchnic.Δσ => 0.0,
  EVtone_Splanchnic.d.x => reflex_delay_init,
  EVtone_Leg.Δσ => 0.0,
  EVtone_Leg.d.x => reflex_delay_init,

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

display(plot(Sol, idxs=[CentralResp.y_A, CentralResp.y_f, PeripheralResp.y_A, PeripheralResp.y_f],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))
display(plot(Sol, idxs=[IschVeins.ΔΘO₂ₛⱼ, IschArterioles.ΔΘO₂ₛⱼ],
        label = ["IschVeins" "IschArterioles" "IschHeart"],
        xlabel = "Time (s)",
        ylabel = "Ischemic Response",
        title = "CNS Ischemic Response"))

display(plot(Sol, idxs=[Abd_A.C.out.cCO₂, Renal_art.C.out.cCO₂, Splanchnic_art.C.out.cCO₂, Leg_art.C.out.cCO₂],
        label = ["Abdomen" "Renal" "Splanchnic" "Leg"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[(Wbₛₚ * Efferent.fab) + (Wcₛₚ * Efferent.fapc) + (Wpₛₚ * Efferent.fasr) + (Wrₛₚ * Efferent.fcpr) - Efferent.θₛₕ]))

display(plot(Sol, idxs=[((Wbₛₚ * Efferent.fab) + (Wcₛₚ * Efferent.fapc) + (Wpₛₚ * Efferent.fasr) + (Wrₛₚ * Efferent.fcpr) - Efferent.θₛₕ),((Wbₛₚ * Efferent.fab) + (Wcₛₚ * Efferent.fapc) + (Wrₛₚ * Efferent.fcpr) - Efferent.θₛₕ),Efferent.fₛₚ],
        label = ["CPR" "No CPR" "Output"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[LungStretchReceptors.fasr]))

display(plot(Sol, idxs=[(Wbₛₚ * Efferent.fab), (Wcₛₚ * Efferent.fapc), (Wpₛₚ * Efferent.fasr), (Wrₛₚ * Efferent.fcpr)]))

display(plot(Sol, idxs=[ABR.fab, Efferent.fₛₚ, EResistance_Leg.Δσ],
        label = ["ABR" "Efferent" "EResistance"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[Cor_art.C.MO₂dyn]))

display(plot(Sol, idxs=[EResistance_Leg.σθ, EResistance_Leg.u],
        label = ["EResistance σθ" "EResistance Δσ"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[EResistance_UpBd.Δσ, EResistance_Renal.Δσ, EResistance_Splanchnic.Δσ, EResistance_Leg.Δσ],
        label = ["EResistance" "EVtone_UpBd" "EVtone_Renal" "EVtone_Splanchnic" "EVtone_Leg"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[RespMuscles.Δp, RespMuscles.BreathInt_held],
        label = ["Peripheral" "Central"],
        xlabel = "Time (s)",
        ylabel = "Efferent",
        title = "Efferent Pathways"))

display(plot(Sol, idxs=[ELH.Δσ, ERH.Δσ],
        label = ["ELH" "ERH"],
        xlabel = "Time (s)",
        ylabel = "Sino-Atrial Node",
        title = "Sino-Atrial Node"))

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
             title = "Arterial Pressure")
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

plot(beat_times, [(Head_art_Vmean + Head_veins_Vmean + Jugular_vein_Vmean + CommonCarotid_Vmean),
        (UpBd_art_Vmean + UpBd_vein_Vmean),
        (Renal_art_Vmean + Renal_vein_Vmean),
        (Splanchnic_art_Vmean + Splanchnic_vein_Vmean),
        (Leg_art_Vmean + Leg_vein_Vmean),
        (Cor_art_Vmean + Cor_vein_Vmean),
        (Asc_A_Vmean + BC_A_Vmean + Thor_A_Vmean + Abd_A_Vmean + Abd_veins_Vmean + Thor_IVC_Vmean + SVC_Vmean),
        (RA_Vmean + RV_Vmean + Pulm_art_Vmean + Pulm_vein_Vmean + LA_Vmean + LV_Vmean)],
        label = ["Head" "Upper Body" "Renal" "Splanchnic" "Leg" "Coronary" "Thoracic" "Cardiopulmonary"],
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Average Branch Volumes")

#### Pulmonary and Respiratory Plots

p4a = plot(Sol, idxs=[RespMuscles.out.p], xlims = (0, 250),
        label = "pₘᵤₛ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Respiratory Muscle Pressure")

p4b = plot(Sol, idxs=[Lungs.pₚₗ], xlims = (0, 250),
        label = "pₚₗ",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Pleural Pressure")

p4c = plot(Sol, idxs=[Lungs.p_A], xlims = (0, 250),
        label = "p_A",
        xlabel = "Time (s)",
        ylabel = "Pressure (mmHg)",
        title = "Alveolar Pressure")

p4d = plot(Sol, idxs=[Lungs.Vrᵢₙ], xlims = (0, 250),
        label = "Air Flow",
        xlabel = "Time (s)",
        ylabel = "Flow (ml/s)",
        title = "Air Flow")

p4e = plot(Sol, idxs=[Lungs.V_A + Lungs.V_D], xlims = (0, 250),
        label = "V_L",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Lung Volume")

p4f = plot(Sol, idxs=[Lungs.V_A], xlims = (0, 250),
        label = "V_A",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Alveolar Volume")

p4g = plot(Sol, idxs=[Lungs.V_D], xlims = (0, 250),
        label = "V_D",
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Dead Space Volume")

p4h = plot(Sol, idxs=[Intrathoracic.pth.p], xlims = (0, 250),
      label = "pₜₕ",
      xlabel = "Time (s)",
      ylabel = "Pressure (mmHg)",
      title = "Intrathoracic Pressure")

      # plot(Sol , idxs=[Lungs.VT])

display(plot(p4a,p4b,p4c,p4d,p4e,p4f,p4g,p4h, layout=(4,2), size=(900,600), suptitle="Lungs"))

display(plot(Sol, idxs=[SA.Eabr_rv_held, SA.Eabr_lv_held],
        xlabel = "Time (s)",
        ylabel = "Volume (ml)",
        title = "Lung Volume"))


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