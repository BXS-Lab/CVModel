"""
Some time constants for later use
"""

T = RRₙₒₘ
Tsys = τᵥₛ * sqrt(T)
Tdias = T - Tsys
Tatria = τₐᵥ * sqrt(T)
Tes = τₐₛ * sqrt(T) # Time varying end-systolic time (s)



Era0 = Ed_ra + (Ees_ra - Ed_ra) *( (Tatria <= Tes) * (1 - cos((Tatria / Tes) * pi)) / 2 + (Tatria > Tes) * (Tatria <= (1.5 * Tes)) * (1 + cos(((Tatria - Tes) / Tes) * 2 * pi)) / 2 + (Tatria > (1.5 * Tes)) * 0)

Ela0 = Ed_la + (Ees_la - Ed_la) *( (Tatria <= Tes) * (1 - cos((Tatria / Tes) * pi)) / 2 + (Tatria > Tes) * (Tatria <= (1.5 * Tes)) * (1 + cos(((Tatria - Tes) / Tes) * 2 * pi)) / 2 + (Tatria > (1.5 * Tes)) * 0)

"""
Find tilt angle at t=0
Here we create a dummy angle driver to extract the initial value of the tilt angle.
"""

@named alphaIC = Alpha()
αeqs = equations(alphaIC)
α_eq = only(αeqs)
αexpr = α_eq.rhs
αdefaults = ModelingToolkit.defaults(alphaIC)
αdefaults[t] = 0.0
αval_evalulated = Symbolics.substitute(αexpr, αdefaults)
α_val_numeric = Symbolics.value(αval_evalulated)

"""
Find gravity at t=0
Here we create a dummy gravity driver to extract the initial value of the gravity.
"""

@named gravityIC = Gravity()
gravityeqs = equations(gravityIC)
gravity_eq = only(gravityeqs)
gravityexpr = gravity_eq.rhs
gravitydefaults = ModelingToolkit.defaults(gravityIC)
gravitydefaults[t] = 0.0
gravityval_evalulated = Symbolics.substitute(gravityexpr, gravitydefaults)
gravity_val_numeric = Symbolics.value(gravityval_evalulated)

"""
Find LBNP at t=0
Here we create a dummy LBNP driver to extract the initial value of the LBNP.
"""

@named lbnpIC = LBNP()
lbnpeqs = equations(lbnpIC)
lbnp_eq = only(lbnpeqs)
lbnpexpr = lbnp_eq.rhs
lbnpdefaults = ModelingToolkit.defaults(lbnpIC)
lbnpdefaults[t] = 0.0
lbnpval_evalulated = Symbolics.substitute(lbnpexpr, lbnpdefaults)
lbnp_val_numeric = Symbolics.value(lbnpval_evalulated)

"""
Determine the initial interstitial volume
Here we use the previously extracted initial values to determine the steady state interstitial volume.
"""

VintIC = Vmax_tilt * sin(α_val_numeric) / sin(85 / 180 * π) * (gravity_val_numeric / 9.81) + Vmax_lbnp * -1 * lbnp_val_numeric / 70

"""
Determine the initial tissue pressures
Here we use the previously extracted initial angle and gravity to determine the initial segment tissue pressures. These are used to find the initial compartment pressures.
"""

pt0_UB = ρ_fft * gravity_val_numeric * (rad_UB/100) * cos(α_val_numeric) * Pa2mmHg
pt0_Thor = ρ_fft * gravity_val_numeric * (rad_Thor/100) * cos(α_val_numeric) * Pa2mmHg
pt0_Abd = ρ_fft * gravity_val_numeric * (rad_Abd/100) * cos(α_val_numeric) * Pa2mmHg
pt0_Leg = ρ_fft * gravity_val_numeric * (rad_Leg/100) * cos(α_val_numeric) * Pa2mmHg

"""
Determine the initial hydrostatic pressures
"""

ph0_Asc_A = ρ_b * gravity_val_numeric * (h_Asc_A/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_BC_A = ρ_b * gravity_val_numeric * (h_BC_A/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_UpBd_art = ρ_b * gravity_val_numeric * (h_UpBd_art/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_UpBd_vein = ρ_b * gravity_val_numeric * (h_UpBd_vein/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_SVC = ρ_b * gravity_val_numeric * (h_SVC/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Thor_A = ρ_b * gravity_val_numeric * (h_Thor_A/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Abd_A = ρ_b * gravity_val_numeric * (h_Abd_A/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Renal_art = ρ_b * gravity_val_numeric * (h_Renal_art/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Renal_vein = ρ_b * gravity_val_numeric * (h_Renal_vein/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Splanchnic_art = ρ_b * gravity_val_numeric * (h_Splanchnic_art/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Splanchnic_vein = ρ_b * gravity_val_numeric * (h_Splanchnic_vein/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Leg_art = ρ_b * gravity_val_numeric * (h_Leg_art/100/con_Leg_art) * sin(α_val_numeric) * Pa2mmHg
ph0_Leg_vein = ρ_b * gravity_val_numeric * (h_Leg_vein/100/con_Leg_vein) * sin(α_val_numeric) * Pa2mmHg
ph0_Abd_veins = ρ_b * gravity_val_numeric * (h_Abd_veins/100/con_default) * sin(α_val_numeric) * Pa2mmHg
ph0_Thor_IVC = ρ_b * gravity_val_numeric * (h_Thor_IVC/100/con_default) * sin(α_val_numeric) * Pa2mmHg

"""
Create Pressure Vector
"""

const N_STATE = 23
x = zeros(N_STATE)
x[1] = 90 # Asc_A.C.p
x[2] = 90 # BC_A.C.p
x[3] = 90 # UpBd_art.C.p
x[4] = 5 # UpBd_vein.C.p
x[5] = 5 # SVC.C.p
x[6] = 90 # Thor_A.C.p
x[7] = 90 # Abd_A.C.p
x[8] = 90 # Renal_art.C.p
x[9] = 5 # Renal_vein.C.p
x[10] = 90 # Splanchnic_art.C.p
x[11] = 5 # Splanchnic_vein.C.p
x[12] = 90 # Leg_art.C.p
x[13] = 5 # Leg_vein.C.p
x[14] = 5 # Abd_veins.C.p
x[15] = 90 # Thor_IVC.C.p
x[16] = 5 # RA.p
x[17] = 5 # RV.p (Diastole)
x[18] = 50 # RV.p (End-Systole)
x[19] = 50 # Pulm_art.C.p
x[20] = 5 # Pulm_vein.C.p
x[21] = 5 # LA.p
x[22] = 5 # LV.p (Diastole)
x[23] = 120 # LV.p (End-Systole)

"""
Create External Pressure Vector
"""

p_rel = zeros(N_STATE)
p_rel[1] = pₜₕ + pt0_Thor
p_rel[2] = pₜₕ + pt0_Thor
p_rel[3] = pt0_UB
p_rel[4] = pt0_UB
p_rel[5] = pₜₕ + pt0_Thor
p_rel[6] = pₜₕ + pt0_Thor
p_rel[7] = pt0_Abd
p_rel[8] = pt0_Abd
p_rel[9] = pt0_Abd
p_rel[10] = pt0_Abd
p_rel[11] = pt0_Abd
p_rel[12] = pt0_Leg
p_rel[13] = pt0_Leg
p_rel[14] = pt0_Abd
p_rel[15] = pₜₕ + pt0_Thor
p_rel[16] = pₜₕ #+ pt0_Thor
p_rel[17] = pₜₕ #+ pt0_Thor
p_rel[18] = pₜₕ #+ pt0_Thor
p_rel[19] = pₜₕ #+ pt0_Thor
p_rel[20] = pₜₕ #+ pt0_Thor
p_rel[21] = pₜₕ #+ pt0_Thor
p_rel[22] = pₜₕ #+ pt0_Thor
p_rel[23] = pₜₕ #+ pt0_Thor

"""
Determine the nonlinear effective compliance equations
"""

consp = (π*C_Splanchnic_vein)/(2*vM_Splanchnic_vein)
conll = (π*C_Leg_vein)/(2*vM_Leg_vein)
conab = (π*C_Abd_veins)/(2*vM_Abd_vein)

      function residuals!(F, x)
        # external pressure vector (precomputed)
        p = p_rel

        # Equation 1: Ventricular volume match
        F[1] =((x[22]-p[22])/Ed_lv - (x[23]-p[23])/Ees_lv) - ((x[17]-p[17])/Ed_rv - (x[18]-p[18])/Ees_rv)
        # Equations 2–22: Flows across resistances equal stroke volume
        SV = (x[22]-p[22])/Ed_lv - (x[23]-p[23])/Ees_lv
        F[2]  = SV - Tsys * (x[23] - x[1]) / R_Asc_A
        F[3]  = SV - T * (x[2] - x[1]) / R_BC_A
        F[4]  = (T * (x[3] - x[4]) / R_UpBd_cap) - T * (x[3] - x[2]) / R_UpBd_art
        # F[5]  = SV - T * (x[3] - x[4]) / R_UpBd_cap
        F[5]  = SV - T * (x[4] - x[5]) / R_UpBd_vein
        F[6]  = SV - T * (x[5] - x[16]) / R_SVC
        F[7]  = SV - T * (x[6] - x[1]) / R_Thor_A
        F[8]  = SV - T * (x[7] - x[6]) / R_Abd_A
        F[9] = (T * (x[8] - x[9]) / R_Renal_cap) - T * (x[8] - x[7]) / R_Renal_art
        # F[11] = SV - T * (x[8] - x[9]) / R_Renal_cap
        F[10] = SV - T * (x[9] - x[14]) / R_Renal_vein
        F[11] = (T * (x[10] - x[11]) / R_Splanchnic_cap) - T * (x[10] - x[7]) / R_Splanchnic_art
        # F[14] = SV - T * (x[10] - x[11]) / R_Splanchnic_cap
        F[12] = SV - T * (x[11] - x[14]) / R_Splanchnic_vein
        F[13] = SV - T * (x[12] - x[7]) / R_Leg_art
        F[14] = SV - T * (x[12] - x[13]) / R_Leg_cap
        F[15] = SV - T * (x[13] - x[14]) / R_Leg_vein
        F[16] = SV - T * (x[14] - x[15]) / R_Abd_veins
        F[17] = SV - T * (x[15] - x[16]) / R_Thor_IVC
        F[18] = SV - Tdias * (x[16] - x[17]) / R_tv
        F[19] = SV - Tsys * (x[18] - x[19]) / Rpa
        F[20] = SV - T * (x[19] - x[20]) / Rpc
        F[21] = SV - T * (x[20] - x[21]) / Rpv
        F[22] = SV - Tdias * (x[21] - x[22]) / R_mv
        F[23] = TBV - (
            (x[1]-p[1])*C_Asc_A + v0_Asc_A +
            (x[2]-p[2])*C_BC_A + v0_BC_A +
            (x[3]-p[3])*C_UpBd_art + v0_UpBd_art +
            (x[4]-p[4])*C_UpBd_vein + v0_UpBd_vein +
            (x[5]-p[5])*C_SVC + v0_SVC +
            (x[6]-p[6])*C_Thor_A + v0_Thor_A +
            (x[7]-p[7])*C_Abd_A + v0_Abd_A +
            (x[8]-p[8])*C_Renal_art + v0_Renal_art +
            (x[9]-p[9])*C_Renal_vein + v0_Renal_vein +
            (x[10]-p[10])*C_Splanchnic_art + v0_Splanchnic_art +
            2 * vM_Splanchnic_vein * atan(consp*(x[11]-p[11])) / π + v0_Splanchnic_vein +
            (x[12]-p[12])*C_Leg_art + v0_Leg_art +
            2 * vM_Leg_vein * atan(conll*(x[13]-p[13])) / π + v0_Leg_vein +
            2 * vM_Abd_vein * atan(conab*(x[14]-p[14])) / π + v0_Abd_veins +
            (x[15]-p[15])*C_Thor_IVC + v0_Thor_IVC +
            (x[16]-p[16])/Era0 + v0_ra +
            (x[17]-p[17])/Ed_rv + v0_rv +
            (x[19]-p[19])*Cpa + v0pa +
            (x[20]-p[20])*Cpv + v0pv +
            (x[21]-p[21])/Ela0 + v0_la +
            (x[22]-p[22])/Ed_lv + v0_lv +
            VintIC
        )
    end

using NLsolve
    sol = nlsolve(residuals!, x)
    x_sol = sol.zero
x = x_sol

    TBV - (
            (x_sol[1]-p_rel[1])*C_Asc_A + v0_Asc_A +
            (x_sol[2]-p_rel[2])*C_BC_A + v0_BC_A +
            (x_sol[3]-p_rel[3])*C_UpBd_art + v0_UpBd_art +
            (x_sol[4]-p_rel[4])*C_UpBd_vein + v0_UpBd_vein +
            (x_sol[5]-p_rel[5])*C_SVC + v0_SVC +
            (x_sol[6]-p_rel[6])*C_Thor_A + v0_Thor_A +
            (x_sol[7]-p_rel[7])*C_Abd_A + v0_Abd_A +
            (x_sol[8]-p_rel[8])*C_Renal_art + v0_Renal_art +
            (x_sol[9]-p_rel[9])*C_Renal_vein + v0_Renal_vein +
            (x_sol[10]-p_rel[10])*C_Splanchnic_art + v0_Splanchnic_art +
            2 * vM_Splanchnic_vein * atan(consp*(x_sol[11]-p_rel[11])) / π + v0_Splanchnic_vein +
            (x_sol[12]-p_rel[12])*C_Leg_art + v0_Leg_art +
            2 * vM_Leg_vein * atan(conll*(x_sol[13]-p_rel[13])) / π + v0_Leg_vein +
            2 * vM_Abd_vein * atan(conab*(x_sol[14]-p_rel[14])) / π + v0_Abd_veins +
            (x_sol[15]-p_rel[15])*C_Thor_IVC + v0_Thor_IVC +
            (x_sol[16]-p_rel[16])/Era0 + v0_ra +
            (x_sol[17]-p_rel[17])/Ed_rv + v0_rv +
            (x_sol[19]-p_rel[19])*Cpa + v0pa +
            (x_sol[20]-p_rel[20])*Cpv + v0pv +
            (x_sol[21]-p_rel[21])/Ela0 + v0_la +
            (x_sol[22]-p_rel[22])/Ed_lv + v0_lv +
            VintIC)

            SV = (x[22]-p_rel[22])/Ed_lv - (x[23]-p_rel[23])/Ees_lv
#     typeof.((Ed_lv, p_rel[22], x[22]))  # Check types







# @parameters p_rel[1:23]
# eqs = [
#     Eq1,
#     Eq2,
#     Eq3,
#     Eq4,
#     Eq5,
#     Eq6,
#     Eq7,
#     Eq8,
#     Eq9,
#     Eq10,
#     Eq11,
#     Eq12,
#     Eq13,
#     Eq14,
#     Eq15,
#     Eq16,
#     Eq17,
#     Eq18,
#     Eq19,
#     Eq20,
#     Eq21,
#     Eq22a,
#     Eq22b,
#     Eq22c,
#     Eq23
# ]

# @named nlsys = NonlinearSystem(eqs, x, p_rel)

# simplified = structural_simplify(nlsys)



# """
# test
# """




# # Constraint Equation:
# tbv_eq = TBVres ~ TBV - (
#     V_ra + V_rv + V_la + V_lv +
#     V_pa + V_pv +
#     V_asc_a + V_bc_a + V_upbd_a + V_thor_a + V_abd_a + V_renal_a + V_spln_a + V_leg_a +
#     V_upbd_v + V_svc + V_renal_v + V_spln_v + V_leg_v + V_abd_v + V_thor_ivc
# )
















# simplified_circ_eqs = [
#   #### Heart and Pulmonary System
#   connect(RA.out, R_tricuspid.in),
#   connect(R_tricuspid.out, RV.in),
#   connect(RV.out, Pulm_art.in),
#   connect(Pulm_art.out, Pulm_cap.in),
#   connect(Pulm_cap.out, Pulm_vein.in),
#   connect(Pulm_vein.out, LA.in),
#   connect(LA.out, R_mitral.in),
#   connect(R_mitral.out, LV.in),
#   connect(LV.out, Asc_A.in),

#   #### Arterial Tree
#   connect(Asc_A.out, BC_A.in, Thor_A.in),
#   connect(BC_A.out, UpBd_art.in),
#   connect(Thor_A.out, Abd_A.in),
#   connect(Abd_A.out, Renal_art.in, Splanchnic_art.in, Leg_art.in),

#   #### Upper Body Circulation
#   connect(UpBd_art.out, UpBd_cap.in),
#   connect(UpBd_cap.out, UpBd_vein.in),

#   #### Renal Circulation
#   connect(Renal_art.out, Renal_cap.in),
#   connect(Renal_cap.out, Renal_vein.in),

#   #### Splanchnic Circulation
#   connect(Splanchnic_art.out, Splanchnic_cap.in),
#   connect(Splanchnic_cap.out, Splanchnic_vein.in),

#   #### Leg Circulation
#   connect(Leg_art.out, Leg_cap.in),
#   connect(Leg_cap.out, Leg_vein.in),

#   #### Venous Tree
#   connect(UpBd_vein.out, SVC.in),
#   connect(Abd_veins.in, Renal_vein.out, Splanchnic_vein.out, Leg_vein.out),
#   connect(Abd_veins.out, Thor_IVC.in),
#   connect(Thor_IVC.out, SVC.out, RA.in),

#   #### Heart Tissue Pressures
#   connect(RA.ep, RA_tissue.out),
#   connect(RV.ep, RV_tissue.out),
#   connect(LA.ep, LA_tissue.out),
#   connect(LV.ep, LV_tissue.out),

#   #### External Pressures
#   connect(Intrathoracic.pth, Asc_A.ep, BC_A.ep, Thor_A.ep, SVC.ep, Thor_IVC.ep, RA_tissue.in, RV_tissue.in, Pulm_art.ep, Pulm_vein.ep, LA_tissue.in, LV_tissue.in),
#   connect(Abdominal.pabd, Abd_A.ep, Renal_art.ep, Splanchnic_art.ep, Renal_vein.ep, Splanchnic_vein.ep, Abd_veins.ep),
#   connect(External.pext, UpBd_art.ep, UpBd_vein.ep),
#   connect(ExternalLBNP.pext, Leg_art.ep, Leg_vein.ep),

#   #### Interstitial Connections (Direct Connections)
#   Splanchnic_vein.C.qint ~ 0.0,
#   Leg_vein.C.qint ~ 0.0,
#   Abd_veins.C.qint ~ 0.0,

#   #### Tilt Equations (Direct Connections)
#   Asc_A.α ~ α_val_numeric,
#   BC_A.α ~ α_val_numeric,
#   UpBd_art.α ~ α_val_numeric,
#   Thor_A.α ~ α_val_numeric,
#   Abd_A.α ~ α_val_numeric,
#   Renal_art.α ~ α_val_numeric,
#   Splanchnic_art.α ~ α_val_numeric,
#   Leg_art.α ~ α_val_numeric,
#   UpBd_vein.α ~ α_val_numeric,
#   SVC.α ~ α_val_numeric,
#   Renal_vein.α ~ α_val_numeric,
#   Splanchnic_vein.α ~ α_val_numeric,
#   Leg_vein.α ~ α_val_numeric,
#   Abd_veins.α ~ α_val_numeric,
#   Thor_IVC.α ~ α_val_numeric,
#   Interstitial.α ~ α_val_numeric,
#   Intrathoracic.α ~ α_val_numeric,
#   Pulm_art.α ~ α_val_numeric,
#   Pulm_vein.α ~ α_val_numeric,
#   RA_tissue.α ~ α_val_numeric,
#   RV_tissue.α ~ α_val_numeric,
#   LA_tissue.α ~ α_val_numeric,
#   LV_tissue.α ~ α_val_numeric,
#   Pulm_cap.α ~ α_val_numeric,

#   #### Gravity Equations (Direct Connections)
#   Asc_A.g ~ gravity_val_numeric,
#   BC_A.g ~ gravity_val_numeric,
#   UpBd_art.g ~ gravity_val_numeric,
#   Thor_A.g ~ gravity_val_numeric,
#   Abd_A.g ~ gravity_val_numeric,
#   Renal_art.g ~ gravity_val_numeric,
#   Splanchnic_art.g ~ gravity_val_numeric,
#   Leg_art.g ~ gravity_val_numeric,
#   UpBd_vein.g ~ gravity_val_numeric,
#   SVC.g ~ gravity_val_numeric,
#   Renal_vein.g ~ gravity_val_numeric,
#   Splanchnic_vein.g ~ gravity_val_numeric,
#   Leg_vein.g ~ gravity_val_numeric,
#   Abd_veins.g ~ gravity_val_numeric,
#   Thor_IVC.g ~ gravity_val_numeric,
#   Interstitial.g ~ gravity_val_numeric,
#   Intrathoracic.g ~ gravity_val_numeric,
#   Pulm_art.g ~ gravity_val_numeric,
#   Pulm_vein.g ~ gravity_val_numeric,
#   RA_tissue.g ~ gravity_val_numeric,
#   RV_tissue.g ~ gravity_val_numeric,
#   LA_tissue.g ~ gravity_val_numeric,
#   LV_tissue.g ~ gravity_val_numeric,
#   Pulm_cap.g ~ gravity_val_numeric,

#   #### LBNP Equations (Direct Connections)
#   ExternalLBNP.p_lbnp ~ lbnp_val_numeric,
#   Interstitial.p_lbnp ~ lbnp_val_numeric,

#   #### Reflex Arc Efferents
#   # The outputs of the transfer functions are connected directly to the relevant compartments via static gains defined in the parameters file.

#   # Arterial Resistance (ABR and CPR)
#   UpBd_cap.R ~ R_UpBd_cap,
#   Renal_cap.R ~ R_Renal_cap,
#   Splanchnic_cap.R ~ R_Splanchnic_cap,
#   Leg_cap.R ~ R_Leg_cap,

#   # Venous Tone (ABR and CPR)
#   UpBd_vein.C.Vabr ~ 0.0,
#   Renal_vein.C.Vabr ~ 0.0,
#   Splanchnic_vein.C.Vabr ~ 0.0,
#   Leg_vein.C.Vabr~ 0.0,
#   UpBd_vein.C.Vcpr ~ 0.0,
#   Renal_vein.C.Vcpr ~ 0.0,
#   Splanchnic_vein.C.Vcpr ~ 0.0,
#   Leg_vein.C.Vcpr ~ 0.0,

#   # Ventricular Contractility (ABR)
#   SA.Eabr_rv ~ 0.0,
#   SA.Eabr_lv ~ 0.0,

#   # Heart Rate (ABR)
#   SA.RRabr ~ 0.0,

#   #### Sino-Atrial Node Connections
#   # The SA node holds the RR Interval and Contractility adjustments, only updating at the start of a new cardiac cycle. The φ signal is the modulated cardiac cycle (atria are offset to contract before the ventricles). The τ signal is the beat-held instantaneous RR interval.
#   RV.Eabr_held ~ 0.0,
#   LV.Eabr_held ~ 0.0,
#   RA.ϕ ~ τₐᵥ * sqrt(RRₙₒₘ),
#   RV.ϕ ~ 0,
#   LA.ϕ ~ τₐᵥ * sqrt(RRₙₒₘ),
#   LV.ϕ ~ 0,
#   RA.τ ~ RRₙₒₘ,
#   RV.τ ~ RRₙₒₘ,
#   LA.τ ~ RRₙₒₘ,
#   LV.τ ~ RRₙₒₘ
# ]

# # Create the simplified circulation model with the volume equation included
# @named _simp_circ_model = ODESystem([simplified_circ_eqs; tbv_eq], t)

# # Compose the main model with all components
# @named simp_circ_model = compose(_simp_circ_model, [RA, R_tricuspid, RV, Pulm_art, Pulm_cap, Pulm_vein, LA, R_mitral, LV, # Heart and Lungs
#   Asc_A, BC_A, UpBd_art, Thor_A, Abd_A, Renal_art, Splanchnic_art, Leg_art, # Arterial Tree
#   UpBd_vein, SVC, Renal_vein, Splanchnic_vein, Leg_vein, Abd_veins, Thor_IVC, # Venous Tree
#   UpBd_cap, Renal_cap, Splanchnic_cap, Leg_cap, # Microcirculation
#   Intrathoracic, Abdominal, External, ExternalLBNP, # External Pressures
#   RA_tissue, RV_tissue, LA_tissue, LV_tissue, # Heart Tissue Pressures
#   alpha_driver, gravity_driver, lbnp_driver # Design of Experiments Drivers
#   ])

# # Now structurally simplify the model
# simp_circ_sys = structural_simplify(simp_circ_model)

# equations(simp_circ_sys)

# param_defaults = Dict()

# for p in parameters(simp_circ_model)
#     pname = ModelingToolkit.getname(p)
#     if hasproperty(ModelParams, pname)
#         param_defaults[p] = getproperty(ModelParams, pname)
#     else
#         # Skip parameters not found in ModelParams
#         continue
#     end
# end

# param_defaults

# param_syms = parameters(simp_circ_model)
# subs_param = Dict()

# for p in param_syms
#     # Convert parameter like :RA₊E to symbol RA₊E and try to look it up in ModelParams
#     pname = nameof(p)
#     if hasproperty(ModelParams, pname)
#         subs_param[p] = getproperty(ModelParams, pname)
#     else
#         @warn "Parameter not found in ModelParams: $pname"
#     end
# end








# # Create a new variable to balance the system
# @variables TBV_param(t)
# tbv_param_eq = TBV_param ~ TBV
# # Create a system with both the volume constraint and parameter definition
# @named volume_constraint = ODESystem([tbv_eq, tbv_param_eq], t)
# @named full_sys = compose(simp_circ_sys, volume_constraint,
#                             [unknowns(simp_circ_sys); TBV_param],
#                             parameters(simp_circ_sys);
#                             observed = observed(simp_circ_sys),
#                             name = :full_sys)



# equations(expand(simp_circ_sys))
# unknowns(simp_circ_sys)
# equations(expand(simp_circ_sys))


# unknowns(volume_constraint)


# @named full_sys = compose(simp_circ_sys, volume_constraint)


# @named full_ss = structural_simplify(full_sys)




# # Get all equations
# eqs = equations(full_sys)
# t = ModelingToolkit.get_iv(full_sys)
# D = Differential(t)

# function extract_dvars(eqs)
#   dvars = Symbolics.Num[]
#   for eq in eqs
#       lhs = eq.lhs
#       if typeof(lhs) == SymbolicUtils.BasicSymbolic{Real}
#           if typeof(lhs.f) == ModelingToolkit.Differential && length(lhs.arguments) == 1
#               push!(dvars, lhs.arguments[1])
#           end
#       end
#   end
#   return dvars
# end

# dvars = extract_dvars(eqs)

# # Identify derivative terms and zero them
# deriv_map = Dict(d => 0.0 for d in dvars)

# # Apply substitution
# ss_eqs = substitute.(eqs, Ref(deriv_map))
# ss_eqs_with_tbv = vcat(ss_eqs, [tbv_eq])

# observed(simp_circ_sys)
# # (Optional) create a new ODESystem if needed
# @named ss_sys = ODESystem(ss_eqs_with_tbv, t, unknowns(simp_circ_sys), parameters(simp_circ_sys))

# using ModelingToolkit, NLsolve

# # Extract symbolic unknowns and equations
# @assert issymbolic(ss_sys)
# vars = unknowns(ss_sys)
# eqs = equations(ss_sys)

# f_expr = [eq.lhs - eq.rhs for eq in eqs]
# f_func = Symbolics.build_function(f_expr, vars; expression=Val{false})[1]
# f_compiled = eval(f_func)

# function f!(F, x)
#     F[:] = f_compiled(x)
# end

# for obs in observed(simp_circ_sys)
#   lhs = obs.lhs
#   @eval const $(Symbol(lhs)) = lhs
# end

# # function f!(F, x)
# #   subs = Dict(zip(vars, x))  # Map variable symbols to values
# #   for (i, eq) in enumerate(eqs)
# #     F[i] = Symbolics.evaluate(Symbolics.substitute(eq.lhs - eq.rhs, subs))
# #   end
# # end

# zip(vars, x)

# using NLsolve
# sol = nlsolve(f!, x0)
# @assert sol.success
# sol_map = Dict(zip(vars, sol.zero))






















# obs_eqs = observed(simp_circ_sys)
# subs_map = Dict(eq.lhs => eq.rhs for eq in obs_eqs)

# eqs_subbed = substitute.(equations(ss_sys), Ref(subs_map))






# observed_subs = Dict(
#     Asc_A₊p_rel => pₜₕ + pt0_Thor,
#     BC_A₊p_rel => pₜₕ + pt0_Thor,
#     UpBd_art₊p_rel => pt0_UB,
#     Thor_A₊p_rel => pₜₕ + pt0_Thor,
#     Abd_A₊p_rel => pt0_Abd,
#     Renal_art₊p_rel => pt0_Abd,
#     Splanchnic_art₊p_rel => pt0_Abd,
#     Leg_art₊p_rel => pt0_Leg,
#     UpBd_vein₊p_rel => pt0_UB,
#     SVC₊p_rel => pₜₕ + pt0_Thor,
#     Renal_vein₊p_rel => pt0_Abd,
#     Splanchnic_vein₊p_rel => pt0_Abd,
#     Leg_vein₊p_rel => pt0_Leg,
#     Abd_veins₊p_rel => pt0_Abd,
#     Thor_IVC₊p_rel => pₜₕ + pt0_Thor,
#     RA₊p_rel => pₜₕ + pt0_Thor,
#     RV₊p_rel => pₜₕ + pt0_Thor,
#     Pulm_art₊C₊p_rel => pₜₕ + pt0_Thor,
#     Pulm_vein₊C₊p_rel => pₜₕ + pt0_Thor,
#     LA₊p_rel => pₜₕ + pt0_Thor,
#     LV₊p_rel => pₜₕ + pt0_Thor,

#     Asc_A₊p_relˍt => 0.0,
#     BC_A₊p_relˍt => 0.0,
#     UpBd_art₊p_relˍt => 0.0,
#     Thor_A₊p_relˍt => 0.0,
#     Abd_A₊p_relˍt => 0.0,
#     Renal_art₊p_relˍt => 0.0,
#     Splanchnic_art₊p_relˍt => 0.0,
#     Leg_art₊p_relˍt => 0.0,
#     UpBd_vein₊p_relˍt => 0.0,
#     SVC₊p_relˍt => 0.0,
#     Renal_vein₊p_relˍt => 0.0,
#     Splanchnic_vein₊p_relˍt => 0.0,
#     Leg_vein₊p_relˍt => 0.0,
#     Abd_veins₊p_relˍt => 0.0,
#     Thor_IVC₊p_relˍt => 0.0,
#     RA₊p_relˍt => 0.0,
#     RV₊p_relˍt => 0.0,
#     Pulm_art₊C₊p_relˍt => 0.0,
#     Pulm_vein₊C₊p_relˍt => 0.0,
#     LA₊p_relˍt => 0.0,
#     LV₊p_relˍt => 0.0,

#     Asc_A₊C₊V₀eff => v0_Asc_A,
#     BC_A₊C₊V₀eff => v0_BC_A,
#     UpBd_art₊C₊V₀eff => v0_UpBd_art,
#     Thor_A₊C₊V₀eff => v0_Thor_A,
#     Abd_A₊C₊V₀eff => v0_Abd_A,
#     Renal_art₊C₊V₀eff => v0_Renal_art,
#     Splanchnic_art₊C₊V₀eff => v0_Splanchnic_art,
#     Leg_art₊C₊V₀eff => v0_Leg_art,
#     UpBd_vein₊C₊V₀eff => v0_UpBd_vein,
#     SVC₊C₊V₀eff => v0_SVC,
#     Renal_vein₊C₊V₀eff => v0_Renal_vein,
#     Splanchnic_vein₊C₊V₀eff => v0_Splanchnic_vein,
#     Leg_vein₊C₊V₀eff => v0_Leg_vein,
#     Abd_veins₊C₊V₀eff => v0_Abd_veins,
#     Thor_IVC₊C₊V₀eff => v0_Thor_IVC,
#     RA₊V₀ => v0_ra,
#     RV₊V₀ => v0_rv,
#     Pulm_art₊C₊V₀eff => v0pa,
#     Pulm_vein₊C₊V₀eff => v0pv,
#     LA₊V₀ => v0_la,
#     LV₊V₀ => v0_lv,

#     Asc_A₊C₊V₀eff_t => 0.0,
#     BC_A₊C₊V₀eff_t => 0.0,
#     UpBd_art₊C₊V₀eff_t => 0.0,
#     Thor_A₊C₊V₀eff_t => 0.0,
#     Abd_A₊C₊V₀eff_t => 0.0,
#     Renal_art₊C₊V₀eff_t => 0.0,
#     Splanchnic_art₊C₊V₀eff_t => 0.0,
#     Leg_art₊C₊V₀eff_t => 0.0,
#     UpBd_vein₊C₊V₀eff_t => 0.0,
#     SVC₊C₊V₀eff_t => 0.0,
#     Renal_vein₊C₊V₀eff_t => 0.0,
#     Splanchnic_vein₊C₊V₀eff_t => 0.0,
#     Leg_vein₊C₊V₀eff_t => 0.0,
#     Abd_veins₊C₊V₀eff_t => 0.0,
#     Thor_IVC₊C₊V₀eff_t => 0.0,
#     RA₊V₀_t => 0.0,
#     RV₊V₀_t => 0.0,
#     Pulm_art₊C₊V₀eff_t => 0.0,
#     Pulm_vein₊C₊V₀eff_t => 0.0,
#     LA₊V₀_t => 0.0,
#     LV₊V₀_t => 0.0
# )

# eqs_sub = substitute.(equations(ss_sys), Ref(observed_subs))



# ss_sys

# prob = ODEProblem(ss_sys, 0)

# @time Sol = solve(ss_sys, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

# vars = unknowns(ss_sys)

# eqs = equations(ss_sys)
# t0 = 0.0  # Evaluate at t = 0

# f_expr = [eq.lhs - eq.rhs for eq in eqs]
# f_func = Symbolics.build_function(f_expr, vars; expression=Val{false})[1]
# f_compiled = eval(f_func)

# function f!(F, x)
#     F[:] = f_compiled(x)
# end

# # function f!(F, x)
# #   subs = Dict(zip(vars, x))  # Map variable symbols to values
# #   for (i, eq) in enumerate(eqs)
# #     F[i] = Symbolics.evaluate(Symbolics.substitute(eq.lhs - eq.rhs, subs))
# #   end
# # end

# zip(vars, x)

# using NLsolve
# ic0 = fill(100.0, length(vars))  # Initial guess, adjust if needed
# F = zeros(length(eqs))

# res = nlsolve(f!, ic0)










# solution = Dict(zip(vars, res.zero))

# p_syms = parameters(ss_sys)

# all_para_names = names(ModelParams, all=true) # Extract all parameter names from the ModelParameters.jl file
# params_syms = filter(n -> isdefined(ModelParams, n) && getfield(ModelParams, n) isa Number, all_para_names) # Filter for numerical (scalar) parameters
# parameter_names  = String.(params_syms) # Convert symbols to strings for DataFrame labeling
# parameter_values = [ getfield(ModelParams, s) for s in params_syms ] # Retrieve corresponding parameter values
# df_parameters = DataFrame(Parameter = parameter_names, Value = parameter_values) # Construct a DataFrame of parameter names and values























# circ_model
# circ_sys

# # Get the expanded equations
# expanded_eqs = equations(expand_connections(circ_model))
# expand_connections(circ_model)

# @variables t
# @variables alpha_driver₊α(t) gravity_driver₊g(t) lbnp_driver₊p_lbnp(t)

# # Set initial values for the drivers
# driver_map = Dict(
#   alpha_driver₊α => α_val_numeric,     # Initial alpha value
#   gravity_driver₊g => gravity_val_numeric,  # Initial gravity value
#   lbnp_driver₊p_lbnp => lbnp_val_numeric,  # Initial LBNP pressure

#   # Interstitial Compartment
#   Interstitial.Qint => 0.0,
#   Interstitial.Vint => VintIC,

#   # Reflex Afferents
#   CPRafferent.x => x0,
#   ABRafferent.x => x0,
#     # Reflex Transfer Functions
#     abr_αr.tfdelay.tftime.x => reflex_delay_init,
#     abr_αr.tfdelay.double_integrator.v => 0.0,
#     abr_αr.tfdelay.double_integrator.y => 0.0,
#     abr_αr.tfpeak.tftime.x => reflex_delay_init,
#     abr_αr.tfpeak.double_integrator.v => 0.0,
#     abr_αr.tfpeak.double_integrator.y => 0.0,
#     abr_αr.tfend.tftime.x => reflex_delay_init,
#     abr_αr.tfend.double_integrator.v => 0.0,
#     abr_αr.tfend.double_integrator.y => 0.0,
#     abr_αv.tfdelay.tftime.x => reflex_delay_init,
#     abr_αv.tfdelay.double_integrator.v => 0.0,
#     abr_αv.tfdelay.double_integrator.y => 0.0,
#     abr_αv.tfpeak.tftime.x => reflex_delay_init,
#     abr_αv.tfpeak.double_integrator.v => 0.0,
#     abr_αv.tfpeak.double_integrator.y => 0.0,
#     abr_αv.tfend.tftime.x => reflex_delay_init,
#     abr_αv.tfend.double_integrator.v => 0.0,
#     abr_αv.tfend.double_integrator.y => 0.0,
#     abr_β.tfdelay.tftime.x => reflex_delay_init,
#     abr_β.tfdelay.double_integrator.v => 0.0,
#     abr_β.tfdelay.double_integrator.y => 0.0,
#     abr_β.tfpeak.tftime.x => reflex_delay_init,
#     abr_β.tfpeak.double_integrator.v => 0.0,
#     abr_β.tfpeak.double_integrator.y => 0.0,
#     abr_β.tfend.tftime.x => reflex_delay_init,
#     abr_β.tfend.double_integrator.v => 0.0,
#     abr_β.tfend.double_integrator.y => 0.0,
#     abr_para.tfdelay.tftime.x => reflex_delay_init,
#     abr_para.tfdelay.double_integrator.v => 0.0,
#     abr_para.tfdelay.double_integrator.y => 0.0,
#     abr_para.tfpeak.tftime.x => reflex_delay_init,
#     abr_para.tfpeak.double_integrator.v => 0.0,
#     abr_para.tfpeak.double_integrator.y => 0.0,
#     abr_para.tfend.tftime.x => reflex_delay_init,
#     abr_para.tfend.double_integrator.v => 0.0,
#     abr_para.tfend.double_integrator.y => 0.0,
#     cpr_αr.tfdelay.tftime.x => reflex_delay_init,
#     cpr_αr.tfdelay.double_integrator.v => 0.0,
#     cpr_αr.tfdelay.double_integrator.y => 0.0,
#     cpr_αr.tfpeak.tftime.x => reflex_delay_init,
#     cpr_αr.tfpeak.double_integrator.v => 0.0,
#     cpr_αr.tfpeak.double_integrator.y => 0.0,
#     cpr_αr.tfend.tftime.x => reflex_delay_init,
#     cpr_αr.tfend.double_integrator.v => 0.0,
#     cpr_αr.tfend.double_integrator.y => 0.0,
#     cpr_αv.tfdelay.tftime.x => reflex_delay_init,
#     cpr_αv.tfdelay.double_integrator.v => 0.0,
#     cpr_αv.tfdelay.double_integrator.y => 0.0,
#     cpr_αv.tfpeak.tftime.x => reflex_delay_init,
#     cpr_αv.tfpeak.double_integrator.v => 0.0,
#     cpr_αv.tfpeak.double_integrator.y => 0.0,
#     cpr_αv.tfend.tftime.x => reflex_delay_init,
#     cpr_αv.tfend.double_integrator.v => 0.0,
#     cpr_αv.tfend.double_integrator.y => 0.0,

#     #### Sino-Atrial Node
#     SA.RR_held => RRₙₒₘ,
#     SA.ϕ => 0.0,
#     SA.Eabr_rv_held => 0.0,
#     SA.Eabr_lv_held => 0.0
# )

# # Apply substitution to the filtered equations
# subbed_eqs = [substitute(eq, driver_map) for eq in expanded_eqs]

# simplified_model = ODESystem(subbed_eqs, t, name=:simplified_model)
# structurally_simplified_model = structural_simplify(simplified_model)

# # Function to check if an equation is reflex-related
# function is_reflex_eq(eq)
#   eq_str = string(eq)
#   reflex_terms = ["cpr_", "abr_", "SA₊", "ABRafferent", "CPRafferent", "Interstitial", "reflex",
#                  "α", "para", "β", "tftime", "tf", "double_integrator"]
#   return any(term -> occursin(term, eq_str), reflex_terms)
# end

# # Filter out reflex-related equations
# filtered_eqs = filter(!is_reflex_eq, expanded_eqs)

# # Print remaining equations
# println("Total number of equations: ", length(expanded_eqs))
# println("Number of equations after filtering reflexes: ", length(filtered_eqs))
# println("\nNon-reflex equations from expanded connections:")
# for (i, eq) in enumerate(filtered_eqs)
#   println("[$i] $eq")
# end

# # Create a substitution map for the driver variables
# @variables t
# @variables alpha_driver₊α(t) gravity_driver₊g(t) lbnp_driver₊p_lbnp(t)

# # Set initial values for the drivers
# driver_map = Dict(
#   alpha_driver₊α => 0.0,     # Initial alpha value
#   gravity_driver₊g => 9.81,  # Initial gravity value
#   lbnp_driver₊p_lbnp => 0.0  # Initial LBNP pressure
# )

# # Apply substitution to the filtered equations
# subbed_eqs = [substitute(eq, driver_map) for eq in filtered_eqs]

# # Use structural_simplify on the equations after substitution
# simplified_eqs = [simplify(eq) for eq in subbed_eqs]

# # Further filter to remove trivial equations (where lhs = rhs after substitution)
# non_trivial_eqs = filter(eq -> !isequal(eq.lhs, eq.rhs), simplified_eqs)

# # Create a simplified model for structural simplification
# @variables t
# simplified_model = ODESystem(non_trivial_eqs, t, name=:simplified_model)
# structurally_simplified_model = structural_simplify(simplified_model)
# structurally_simplified_eqs = equations(structurally_simplified_model)

# println("Number of equations before driver substitution: ", length(filtered_eqs))
# println("Number of equations after driver substitution and simplification: ", length(non_trivial_eqs))
# println("Number of equations after structural simplification: ", length(structurally_simplified_eqs))
# println("\nStructurally simplified equations:")
# for (i, eq) in enumerate(structurally_simplified_eqs)
#   println("[$i] $eq")
# end



# parameters(circ_sys)

# expanded = equations(ModelingToolkit.expand_connections(circ_model))
# observed(circ_sys)

# @variables t
# @variables alpha_driver₊α(t) gravity_driver₊g(t) lbnp_driver₊p_lbnp(t)

# # Substitution map for driver variables
# subs_map = Dict(
#   alpha_driver₊α => 0.0,     # Initial alpha value
#   gravity_driver₊g => 9.81,    # Initial gravity value
#   lbnp_driver₊p_lbnp => 0.0  # Initial LBNP pressure
# )

# # Apply substitution map to the expanded equations
# subbed_eqs = [substitute(eq, subs_map) for eq in expanded]

# # Simplify the equations after substitution
# simplified_eqs = [simplify(eq) for eq in subbed_eqs]

# # Let's filter more aggressively to remove alpha-related equations
# filtered_eqs = filter(eq -> begin
#   eq_str = string(eq)

#   # Skip equations with alpha terms since we know alpha_driver₊α = 0
#   if occursin("α", eq_str)
#     return false
#   end

#   # Keep only non-trivial equations (where lhs != rhs)
#   !isequal(eq.lhs, eq.rhs)
# end, simplified_eqs)

# println("Original number of equations: ", length(expanded))
# println("Number of equations after filtering: ", length(filtered_eqs))

# # Use filtered_eqs for subsequent processing instead of simplified_eqs
# non_trivial_eqs = filtered_eqs
# subbed_eqs = [substitute(eq, subs_map) for eq in expanded]

# # Simplify the equations after substitution
# simplified_eqs = [simplify(eq) for eq in subbed_eqs]
# println("Number of equations before filtering: ", length(simplified_eqs))

# # Filter out trivial equations (where lhs equals rhs) but keep α-related equations
# non_trivial_eqs = filter(eq -> begin
#     eq_str = string(eq)
#     # Keep only non-trivial equations that don't contain α
#     is_non_trivial = !isequal(eq.lhs, eq.rhs)
#     has_no_alpha = !occursin("α", eq_str)
#     # Or keep equations with α only if they don't simplify to zero
#     is_alpha_non_zero = occursin("α", eq_str) && !isequal(simplify(eq.lhs - eq.rhs), 0)

#     # Return true to keep the equation
#     (is_non_trivial && has_no_alpha) || is_alpha_non_zero
# end, simplified_eqs)

# println("Number of equations after filtering: ", length(non_trivial_eqs))















# eq_objs = ModelingToolkit.get_eqs(circ_model)
# conns = filter(eq -> eq isa ModelingToolkit.Connection, eq_objs)
# typeof(eq_objs[119])


# conn_eqs = filter(eq -> occursin("connect", string(eq)), ModelingToolkit.get_eqs(circ_model))

# simplified_model = ODESystem(conn_eqs, t, name=:simplified_model)
# structurally_simplified_model = structural_simplify(simplified_model)
# # Print the structurally simplified model equations
# println("Structurally simplified model equations:")
# for (i, eq) in enumerate(equations(structurally_simplified_model))
#   println("[$i] $eq")
# end

# # Print observed variables
# println("\nObserved variables:")
# for (i, var) in enumerate(observed(structurally_simplified_model))
#   println("[$i] $var")
# end







# connections = [eq.rhs for eq in conn_eqs if eq.rhs isa ModelingToolkit.Connection]
# conn_subsystems = unique(Iterators.flatten(c.systems for c in connections))
# @named conn_only_sys = compose(ODESystem[], conn_subsystems...)
# expanded_conn_sys = expand_connections(conn_only_sys)
# eqs = equations(expanded_conn_sys)
# conn_onon_trivial_eqs = filter(eq -> begin
#     eq_str = string(eq)
#     # Keep only non-trivial equations that don't contain α
#     is_non_trivial = !isequal(eq.lhs, eq.rhs)
#     has_no_alpha = !occursin("α", eq_str)p    # Or keep equations with α only if they don't simplify to zero
#     is_alpha_non_zero = occursin("α", eq_str) && !isequal(simplify(eq.lhs - eq.rhs), 0)

#     # Return true to keep the equation
#     (is_non_trivial && has_no_alpha) || is_alpha_non_zero
# end, simplified_eqs)

# println("Number of equations after filtering: ", length(non_trivial_eqs))ha_driver₊α(t) lbnp_driver₊p_lbnp(t)

# ic_map = Dict(
#   gravity_driver₊g(t) => gravity_val_numeric,
#   alpha_driver₊α(t) => 0.0,
#   lbnp_driver₊p_lbnp(t) => lbnp_val_numeric,
#   )  # fill in with your known values

# eqs_sub = substitute.(expanded_eqs, Ref(ic_map))
# eqs_simplified = simplify.(eqs_sub)
# reduced_eqs = filter(eq -> !isequal(eq.lhs, eq.rhs), eqs_simplified)












# for (i, eq) in enumerate(expanded_eqs)
#     println("[$i] ", eq)
# end

# conn = conn_eqs[1].rhs
# subsystems = conn.systems
# vars = reduce(vcat, [
#     filter(v -> ModelingToolkit.getname(v) in (:p, :q), ModelingToolkit.get_variables(sys))
#     for sys in subsystems
# ])












# equations(expand_connections(conn_eqs[1].rhs.systems[1]))

# function extract_connection_ports(eq)
#   rhs = eq.rhs
#   if hasproperty(rhs, :f) && rhs.f == :connect && hasproperty(rhs, :args)
#       return rhs.args  # this is a tuple of connected ports
#   else
#       return nothing
#   end
# end

# function build_flow_pressure_eqs(conn_eqs, t)
#   flow_eqs = Equation[]
#   pressure_eqs = Equation[]

#   for eq in conn_eqs
#       ports = extract_connection_ports(eq)
#       if ports === nothing
#           continue
#       end

#       # Get flow and pressure variables for each port
#       q_vars = [getproperty(p, :q)(t) for p in ports]
#       p_vars = [getproperty(p, :p)(t) for p in ports]

#       # Flow conservation: first port is sum of others
#       push!(flow_eqs, sum(q_vars[2:end]) ~ q_vars[1])

#       # Pressure equality
#       for i in 2:length(p_vars)
#           push!(pressure_eqs, p_vars[i] ~ p_vars[1])
#       end
#   end

#   return flow_eqs, pressure_eqs
# end

# flow_eqs, pressure_eqs = build_flow_pressure_eqs(conn_eqs, t)


# eq = conn_eqs[1]
# typeof(eq)         # should be Equation
# typeof(eq.rhs)     # what's the type of the RHS?
# dump(eq.rhs)

# function extract_connected_vars(conn_eqs)
#   connected_vars = []
#   for eq in conn_eqs
#       if eq.rhs isa ModelingToolkit.Connection
#           conn = eq.rhs
#           # The systems field holds the connected components
#           systems = conn.systems
#           for sys in systems
#               append!(connected_vars, ModelingToolkit.unknowns(sys))
#           end
#       end
#   end
#   return connected_vars
# end

# function build_flow_balance_eqs(conn_eqs, t)
#   flow_eqs = Equation[]
#   for eq in conn_eqs
#       conn = eq.rhs
#       systems = conn.systems
#       inflow = Symbolics.Num[]
#       outflow = Symbolics.Num[]
#       for sys in systems
#           qvars = ModelingToolkit.unknowns(sys)
#           for v in qvars
#               if ModelingToolkit.getname(v) == :q
#                   # Guess direction from connector name
#                   if occursin("in", string(sys.name))
#                       push!(inflow, v)
#                   elseif occursin("out", string(sys.name))
#                       push!(outflow, v)
#                   end
#               end
#           end
#       end
#       if !isempty(inflow) || !isempty(outflow)
#           push!(flow_eqs, sum(outflow) ~ sum(inflow))
#       end
#   end
#   return flow_eqs
# end

# flow_eqs = build_flow_balance_eqs(conn_eqs, t)

# conn_eqs[1].rhs.systems[2]

# unknowns(conn_eqs[1].rhs.systems[2])






























# function flatten_expr(expr)
#   if expr isa Num && expr isa Term && expr.f == +
#       return reduce(vcat, flatten_expr.(expr.args))
#   else
#       return [expr]
#   end
# end

# function is_flow_balance(eq)
#   terms = flatten_expr(eq.lhs - eq.rhs)
#   return any(term -> occursin("_q", string(term)), terms)
# end

# flow_balance_eqs = filter(is_flow_balance, equations(circ_model))
















# function is_flow_balance(eq)
#   terms = Symbolics.ordered_ops(eq.rhs - eq.lhs)
#   return any(x -> occursin("_q", string(x)), terms)
# end

# flow_balance_eqs = filter(is_flow_balance, all_eqs)

# all_eqs = equations(circ_sys)
# connection_eqs = filter(eq -> eq.lhs != eq.rhs && Symbolics.issym(eq.lhs) && Symbolics.issym(eq.rhs), all_eqs)

# function is_connection_eq(eq)
#   lhs = eq.lhs
#   rhs = eq.rhs
#   # Must be equality between two symbols (e.g., q1(t) ~ q2(t))
#   return Symbolics.issym(lhs) && Symbolics.issym(rhs) && lhs != rhs
# end

# connection_eqs = filter(is_connection_eq, equations(circ_model))
# equations(circ_model)











# t = ModelingToolkit.get_iv(circ_sys)
# D = Differential(t)

# # Extract variables of the form D(var) ~ f(...) and return var
# function extract_dvars(eqs)
#   dvars = Symbolics.Num[]
#   for eq in eqs
#       lhs = eq.lhs
#       if typeof(lhs) == SymbolicUtils.BasicSymbolic{Real}
#           if typeof(lhs.f) == ModelingToolkit.Differential && length(lhs.arguments) == 1
#               push!(dvars, lhs.arguments[1])
#           end
#       end
#   end
#   return dvars
# end

# dvars = extract_dvars(eqs)

# deriv_map = Dict(D(var) => 0.0 for var in dvars)
# ss_eqs = map(eq -> substitute(eq, deriv_map), eqs)

# ic_map = Dict(
#    #### Interstitial Compartment
#    Interstitial.Qint => 0.0,
#    Interstitial.Vint => VintIC,

#    #### Reflex Afferents
#    CPRafferent.x => x0,
#    ABRafferent.x => x0,

#    #### Reflex Transfer Functions
#    abr_αr.tfdelay.tftime.x => reflex_delay_init,
#    abr_αr.tfdelay.double_integrator.v => 0.0,
#    abr_αr.tfdelay.double_integrator.y => 0.0,
#    abr_αr.tfpeak.tftime.x => reflex_delay_init,
#    abr_αr.tfpeak.double_integrator.v => 0.0,
#    abr_αr.tfpeak.double_integrator.y => 0.0,
#    abr_αr.tfend.tftime.x => reflex_delay_init,
#    abr_αr.tfend.double_integrator.v => 0.0,
#    abr_αr.tfend.double_integrator.y => 0.0,
#    abr_αv.tfdelay.tftime.x => reflex_delay_init,
#    abr_αv.tfdelay.double_integrator.v => 0.0,
#    abr_αv.tfdelay.double_integrator.y => 0.0,
#    abr_αv.tfpeak.tftime.x => reflex_delay_init,
#    abr_αv.tfpeak.double_integrator.v => 0.0,
#    abr_αv.tfpeak.double_integrator.y => 0.0,
#    abr_αv.tfend.tftime.x => reflex_delay_init,
#    abr_αv.tfend.double_integrator.v => 0.0,
#    abr_αv.tfend.double_integrator.y => 0.0,
#    abr_β.tfdelay.tftime.x => reflex_delay_init,
#    abr_β.tfdelay.double_integrator.v => 0.0,
#    abr_β.tfdelay.double_integrator.y => 0.0,
#    abr_β.tfpeak.tftime.x => reflex_delay_init,
#    abr_β.tfpeak.double_integrator.v => 0.0,
#    abr_β.tfpeak.double_integrator.y => 0.0,
#    abr_β.tfend.tftime.x => reflex_delay_init,
#    abr_β.tfend.double_integrator.v => 0.0,
#    abr_β.tfend.double_integrator.y => 0.0,
#    abr_para.tfdelay.tftime.x => reflex_delay_init,
#    abr_para.tfdelay.double_integrator.v => 0.0,
#    abr_para.tfdelay.double_integrator.y => 0.0,
#    abr_para.tfpeak.tftime.x => reflex_delay_init,
#    abr_para.tfpeak.double_integrator.v => 0.0,
#    abr_para.tfpeak.double_integrator.y => 0.0,
#    abr_para.tfend.tftime.x => reflex_delay_init,
#    abr_para.tfend.double_integrator.v => 0.0,
#    abr_para.tfend.double_integrator.y => 0.0,
#    cpr_αr.tfdelay.tftime.x => reflex_delay_init,
#    cpr_αr.tfdelay.double_integrator.v => 0.0,
#    cpr_αr.tfdelay.double_integrator.y => 0.0,
#    cpr_αr.tfpeak.tftime.x => reflex_delay_init,
#    cpr_αr.tfpeak.double_integrator.v => 0.0,
#    cpr_αr.tfpeak.double_integrator.y => 0.0,
#    cpr_αr.tfend.tftime.x => reflex_delay_init,
#    cpr_αr.tfend.double_integrator.v => 0.0,
#    cpr_αr.tfend.double_integrator.y => 0.0,
#    cpr_αv.tfdelay.tftime.x => reflex_delay_init,
#    cpr_αv.tfdelay.double_integrator.v => 0.0,
#    cpr_αv.tfdelay.double_integrator.y => 0.0,
#    cpr_αv.tfpeak.tftime.x => reflex_delay_init,
#    cpr_αv.tfpeak.double_integrator.v => 0.0,
#    cpr_αv.tfpeak.double_integrator.y => 0.0,
#    cpr_αv.tfend.tftime.x => reflex_delay_init,
#    cpr_αv.tfend.double_integrator.v => 0.0,
#    cpr_αv.tfend.double_integrator.y => 0.0,

#    #### Sino-Atrial Node
#    SA.RR_held => RRₙₒₘ,
#    SA.ϕ => 0.0,
#    SA.Eabr_rv_held => 0.0,
#    SA.Eabr_lv_held => 0.0
# )

# ss_eqs_sub = map(eq -> substitute(eq, ic_map), ss_eqs)

# nontrivial_eqs = filter(eq -> !Symbolics.isequal(Symbolics.simplify(eq.lhs), Symbolics.simplify(eq.rhs)), ss_eqs_sub)

# remaining_eqs = filter(eq -> !issubset(operations(eq), keys(ic_map)), ss_eqs_sub)










# function is_reflex_eq(eq)
#   txt = string(eq)
#   return occursin("cpr_", txt) || occursin("abr_", txt) || occursin("SA₊", txt) || occursin("ABR", txt) || occursin("CPR", txt) || occursin("Interstitial", txt)
# end

# ss_eqs_filtered = filter(!is_reflex_eq, ss_eqs)

# function collect_qint_vars(expr, qint_vars=Symbolics.Num[])
#   if expr isa SymbolicUtils.BasicSymbolic{Real}
#       if occursin("qint", string(expr)) &&
#          hasproperty(expr, :arguments) &&
#          length(expr.arguments) == 1 &&
#          string(expr.arguments[1]) == "t" &&
#          !any(x -> isequal(x, expr), qint_vars)

#           push!(qint_vars, expr)
#       end
#       if hasproperty(expr, :arguments)
#           for subexpr in expr.arguments
#               collect_qint_vars(subexpr, qint_vars)
#           end
#       end
#   end
#   return qint_vars
# end

# function find_qint_vars(eqs)
#   qint_vars = Symbolics.Num[]
#   for eq in eqs
#       collect_qint_vars(eq.rhs, qint_vars)
#   end
#   return qint_vars
# end

# typeof(ss_eqs_filtered[19].rhs)
# qint_vars = find_qint_vars(ss_eqs_filtered)
# qint_map = Dict(q => 0.0 for q in qint_vars)

# ss_eqs_cleaned = map(eq -> substitute(eq, qint_map), ss_eqs_filtered)

# @unpack RA, RV, LA, LV, Pulm_art, Pulm_vein,
# Asc_A, BC_A, UpBd_art, Thor_A, Abd_A, Renal_art, Splanchnic_art, Leg_art,
# UpBd_vein, SVC, Renal_vein, Splanchnic_vein, Leg_vein, Abd_veins, Thor_IVC = circ_sys

# V_ra      = (RA.p - RA.p_rel) / RA.E + RA.V₀
# V_rv      = (RV.p - RV.p_rel) / RV.E + RV.V₀
# V_la      = (LA.p - LA.p_rel) / LA.E + LA.V₀
# V_lv      = (LV.p - LV.p_rel) / LV.E + LV.V₀

# V_pa      = (Pulm_art.C.p - Pulm_art.C.p_rel) * Pulm_art.C.C + Pulm_art.C.V₀eff
# V_pv      = (Pulm_vein.C.p - Pulm_vein.C.p_rel) * Pulm_vein.C.C + Pulm_vein.C.V₀eff

# V_asc_a   = (Asc_A.C.p - Asc_A.C.p_rel) * Asc_A.C.C + Asc_A.C.V₀eff
# V_bc_a    = (BC_A.C.p - BC_A.C.p_rel) * BC_A.C.C + BC_A.C.V₀eff
# V_upbd_a  = (UpBd_art.C.p - UpBd_art.C.p_rel) * UpBd_art.C.C + UpBd_art.C.V₀eff
# V_thor_a  = (Thor_A.C.p - Thor_A.C.p_rel) * Thor_A.C.C + Thor_A.C.V₀eff
# V_abd_a   = (Abd_A.C.p - Abd_A.C.p_rel) * Abd_A.C.C + Abd_A.C.V₀eff
# V_renal_a = (Renal_art.C.p - Renal_art.C.p_rel) * Renal_art.C.C + Renal_art.C.V₀eff
# V_spln_a  = (Splanchnic_art.C.p - Splanchnic_art.C.p_rel) * Splanchnic_art.C.C + Splanchnic_art.C.V₀eff
# V_leg_a   = (Leg_art.C.p - Leg_art.C.p_rel) * Leg_art.C.C + Leg_art.C.V₀eff

# V_upbd_v  = (UpBd_vein.C.p - UpBd_vein.C.p_rel) * UpBd_vein.C.C + UpBd_vein.C.V₀eff
# V_svc     = (SVC.C.p - SVC.C.p_rel) * SVC.C.C + SVC.C.V₀eff
# V_renal_v = (Renal_vein.C.p - Renal_vein.C.p_rel) * Renal_vein.C.C + Renal_vein.C.V₀eff
# V_spln_v  = (Splanchnic_vein.C.p - Splanchnic_vein.C.p_rel) * Splanchnic_vein.C.Cvar + Splanchnic_vein.C.V₀eff
# V_leg_v   = (Leg_vein.C.p - Leg_vein.C.p_rel) * Leg_vein.C.Cvar + Leg_vein.C.V₀eff
# V_abd_v   = (Abd_veins.C.p - Abd_veins.C.p_rel) * Abd_veins.C.Cvar + Abd_veins.C.V₀eff
# V_thor_ivc= (Thor_IVC.C.p - Thor_IVC.C.p_rel) * Thor_IVC.C.C + Thor_IVC.C.V₀eff

# tbv_eq = 0 ~ TBV - (
#     V_ra + V_rv + V_la + V_lv +
#     V_pa + V_pv +
#     V_asc_a + V_bc_a + V_upbd_a + V_thor_a + V_abd_a + V_renal_a + V_spln_a + V_leg_a +
#     V_upbd_v + V_svc + V_renal_v + V_spln_v + V_leg_v + V_abd_v + V_thor_ivc
# )

# push!(ss_eqs_cleaned, tbv_eq)

# all_globals = names(ModelParams, all=true)

# # param_syms = []
# # param_vals = []

# # for name in all_globals
# #   # Skip internals and symbols not in your file
# #   if isdefined(Main, name)
# #       val = getfield(Main, name)
# #       if isa(val, Number)
# #           push!(param_syms, Symbol(name))
# #           push!(param_vals, val)
# #       end
# #   end
# # end

# # params_vector = Pair.(param_syms, param_vals)

# # vars = [
# #     Asc_A.C.p, BC_A.C.p, UpBd_art.C.p, Thor_A.C.p, Abd_A.C.p,
# #     Renal_art.C.p, Splanchnic_art.C.p, Leg_art.C.p,
# #     UpBd_vein.C.p, SVC.C.p, Renal_vein.C.p, Splanchnic_vein.C.p,
# #     Leg_vein.C.p, Abd_veins.C.p, Thor_IVC.C.p,
# #     Pulm_art.C.p, Pulm_vein.C.p,
# #     RA.p, RV.p, LA.p, LV.p
# # ]

# # using NLsolve

# # subs = Dict(params_vector...)
# # eqs_sub = map(eq -> substitute(eq, subs), ss_eqs_cleaned)

# # # Define a function from your equation vector
# # function f!(F, x)
# #     subs = Dict(zip(vars, x))
# #     F[:] = substitute.([eq.rhs - eq.lhs for eq in eqs_sub], subs)
# #   end

# # x0 = ones(length(vars))  # initial guess
# # res = nlsolve(f!, x0)

# # ics = ModelingToolkit.solve_for(ss_eqs_cleaned, vars; parameters=params_vector)









# # t




















# # vars = get_variables(circ_sys)

# # dvars = [v for v in ModelingToolkit.variables(sys) if ModelingToolkit.isdifferential(v)]

# # dvars = get_variables(circ_sys)
# # display(dvars)
# # dmap = Dict(d => 0.0 for d in dvars)

# # ss_eqs = substitute.(eqs, dmap)

# # using .ModelParameters  # to access TBV
# # using .YourModelModules

# # Vtotal_expr = sum(comp.C.V for comp in all_comps_with_volume) + Interstitial.Vint
# # TBV_constraint = Vtotal_expr ~ TBV

# # push!(ss_eqs, TBV_constraint)

# # function collect_volumes(sys)
# #   return [getfield(comp, :C).V for comp in components(sys) if hasfield(typeof(getfield(comp, :C)), :V)]
# # end

# # ss_sys = NonlinearSystem(ss_eqs, ModelingToolkit.states(sys), parameters(sys))

# # using NLsolve
# # u0_guess = [...]  # you must provide a reasonable initial guess
# # p = [...]         # use your model parameters

# # prob = NonlinearProblem(ss_sys, u0_guess, p)
# # sol = solve(prob, NewtonRaphson())

# # sol = ModelingToolkit.solve_for(ss_eqs, ModelingToolkit.states(sys))




# # ss_eqs = [lhs ~ 0 for eq in eqs for lhs in ModelingToolkit.get_differential_vars(eq)]

# # is_derivative_eq(eq) = any(occursin("D(", string(term)) for term in Symbolics.ordered_terms(eq))
# # ss_eqs = [eq => 0.0 for eq in eqs if is_derivative_eq(eq)]

# """
# Determine the constitutive volume equations
# """

# V_asc_a   = (x[1] - p_rel[1]) * C_Asc_A + v0_Asc_A
# V_bc_a    = (x[2] - p_rel[2]) * C_BC_A + v0_BC_A
# V_upbd_a  = (x[3] - p_rel[3]) * C_UpBd_art + v0_UpBd_art
# V_upbd_v  = (x[4] - p_rel[4]) * C_UpBd_vein + v0_UpBd_vein
# V_svc     = (x[5] - p_rel[5]) * C_SVC + v0_SVC
# V_thor_a  = (x[6] - p_rel[6]) * C_Thor_A + v0_Thor_A
# V_abd_a   = (x[7] - p_rel[7]) * C_Abd_A + v0_Abd_A
# V_renal_a = (x[8] - p_rel[8]) * C_Renal_art + v0_Renal_art
# V_renal_v = (x[9] - p_rel[9]) * C_Renal_vein + v0_Renal_vein
# V_spln_a  = (x[10] - p_rel[10]) * Csp + v0_Splanchnic_art
# V_spln_v  = (x[11] - p_rel[11]) * Splanchnic_vein.C.Cvar + v0_Splanchnic_vein
# V_leg_a   = (x[12] - p_rel[12]) * C_Leg_art + v0_Leg_art
# V_leg_v   = (x[13] - p_rel[13]) * Cll + v0_Leg_vein
# V_abd_v   = (x[14] - p_rel[14]) * Cab + v0_Abd_veins
# V_thor_ivc= (x[15] - p_rel[15]) * C_Thor_IVC + v0_Thor_IVC

# V_ra      = (x[16] - p_rel[16]) / Ed_ra + v0_ra
# V_rv      = (x[17] - p_rel[17]) / Ed_rv + v0_rv

# V_pa      = (x[19] - p_rel[19]) * Cpa + v0pa
# V_pv      = (x[20] - p_rel[20]) * Cpv+ v0pv

# V_la      = (x[21] - p_rel[21]) / Ed_la + v0_la
# V_lv      = (x[22] - p_rel[22]) / Ed_lv + v0_lv

# """
# Equation 1:
# The difference between diastolic and systolic volumes are equal in the left and right ventricles
# """

# Eq1 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ (x[17]-p_rel[17]) / Ed_rv - (x[18]-p_rel[18]) / Ees_rv

# """
# Equations 2–22: Equate that volume with the volume which flows through each resistance during the time that it is active
# (Actually 24 equations but simplifies down to 21)
# """

# # R1 (Asc_A)
# Eq2 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ Tsys * (x[1] - x[23]) / R_Asc_A

# # R2 (BC_A)
# Eq3 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[2] - x[1]) / R_BC_A

# # R3 (UpBd_art)
# Eq4 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[3] - x[2]) / R_UpBd_art

# # Rub (UpBd_cap)
# Eq5 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[4] - x[3]) / R_UpBd_cap

# # R4 (UpBd_vein)
# Eq6 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[5] - x[4]) / R_UpBd_vein

# # R5 (SVC)
# Eq7 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[16] - x[5]) / R_SVC

# # R6 (Thor_A)
# Eq8 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[6] - x[1]) / R_Thor_A

# # R7 (Abd_A)
# Eq9 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[7] - x[6]) / R_Abd_A

# # R8 (Renal_art)
# Eq10 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[8] - x[7]) / R_Renal_art

# # Rrc (Renal_cap)
# Eq11 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[9] - x[8]) / R_Renal_cap

# # R9 (Renal_vein)
# Eq12 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[14] - x[9]) / R_Renal_vein

# # R10 (Splanchnic_art)
# Eq13 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[10] - x[7]) / R_Splanchnic_art

# # Rsc (Splanchnic_cap)
# Eq14 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[11] - x[10]) / R_Splanchnic_cap

# # R11 (Splanchnic_vein)
# Eq15 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[14] - x[11]) / R_Splanchnic_vein

# # R12 (Leg_art)
# Eq16 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[12] - x[7]) / R_Leg_art

# # Rlc (Leg_cap)
# Eq17 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[13] - x[12]) / R_Leg_cap

# # R13 (Leg_vein)
# Eq18 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[14] - x[13]) / R_Leg_vein

# # R14 (Abd_veins)
# Eq19 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[15] - x[14]) / R_Abd_veins

# # R15 (Thor_IVC)
# Eq20 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[16] - x[15]) / R_Thor_IVC

# # r_tv (Tricuspid Valve)
# Eq21 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ Tdias * (x[17] - x[16]) / R_tv

# # Rpa (Pulm_art)
# Eq22 = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ Tsys * (x[19] - x[18]) / Rpa

# # Rpc (Pulm_cap)
# Eq22a = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[20] - x[19]) / Rpc

# # Rpv (Pulm_vein)
# Eq22b = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ T * (x[21] - x[20]) / Rpv

# # r_mv (Mitral Valve)
# Eq22c = (x[22]-p_rel[22]) / Ed_lv - (x[23]-p_rel[23]) / Ees_lv ~ Tdias * (x[22] - x[21]) / R_mv

# """
# Equation 23: Total Blood Volume is equal to the sum of all compartment volumes
# """