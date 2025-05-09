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

"""

pₜₕ₀ = pₚₗ₀ - 3.5 * (gravity_val_numeric / 9.81) * sin(α_val_numeric)

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
p_rel[1] = pₜₕ₀ + pt0_Thor
p_rel[2] = pₜₕ₀ + pt0_Thor
p_rel[3] = pt0_UB
p_rel[4] = pt0_UB
p_rel[5] = pₜₕ₀ + pt0_Thor
p_rel[6] = pₜₕ₀ + pt0_Thor
p_rel[7] = pt0_Abd
p_rel[8] = pt0_Abd
p_rel[9] = pt0_Abd
p_rel[10] = pt0_Abd
p_rel[11] = pt0_Abd
p_rel[12] = pt0_Leg
p_rel[13] = pt0_Leg
p_rel[14] = pt0_Abd
p_rel[15] = pₜₕ₀ + pt0_Thor
p_rel[16] = pₜₕ₀
p_rel[17] = pₜₕ₀
p_rel[18] = pₜₕ₀
p_rel[19] = pₜₕ₀
p_rel[20] = pₜₕ₀
p_rel[21] = pₜₕ₀
p_rel[22] = pₜₕ₀
p_rel[23] = pₜₕ₀

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
        F[18] = SV - Tdias * (x[16] - x[17]) / sqrt(1/(ρ_b / (2 * Ann_tv^2)))
        F[19] = SV - Tsys * (x[18] - x[19]) / Rpa
        F[20] = SV - T * (x[19] - x[20]) / Rpc
        F[21] = SV - T * (x[20] - x[21]) / Rpv
        F[22] = SV - Tdias * (x[21] - x[22]) / sqrt(1/(ρ_b / (2 * Ann_mv^2)))
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