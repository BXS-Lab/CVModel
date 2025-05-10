"""
This file determines the initial compartment pressures based on the method developed by Heldt (2004). It is modified to account for external tissue pressures (Whittle, 2023). We also dynamically extract the initial values for the tilt angle, gravity, and LBNP. These are used to set the initial interstitial volume to a steady state value such that Qint = 0.0.
"""

"""
Preamble
This code section sets up the initial variables and extracts the relevant parameters necessary for the routine.
"""

const N_STATE = 23

A = zeros(N_STATE, N_STATE)
b = zeros(N_STATE)

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

pt0_UB = ρ_fft * gravity_val_numeric * rad_UB * cos(α_val_numeric) * 0.0000750062
pt0_Thor = ρ_fft * gravity_val_numeric * rad_Thor * cos(α_val_numeric) * 0.0000750062
pt0_Abd = ρ_fft * gravity_val_numeric * rad_Abd * cos(α_val_numeric) * 0.0000750062
pt0_Leg = ρ_fft * gravity_val_numeric * rad_Leg * cos(α_val_numeric) * 0.0000750062

"""
This function defines a modified Newton Raphson solver to solve the system of equations for the initial compartment pressures. An initial guess is obtained by solving a matrix equation (below) using the zero pressure filling volumes of the nonlinear compartments. The mnewt!() function is subsequently called iteratively to converge on the true pressures when incorporating the nonlinearities. This method was defined by Heldt (2004) and modified by Whittle (2023) to include the external tissue pressures.

Note: Heldt's C code reversed the order of the rows and columns in the A matrix. This is corrected here.
"""

function mnewt!(b)

  consp = (π*C_Splanchnic_vein)/(2*vM_Splanchnic_vein)
  Csp = C_Splanchnic_vein/(1+consp^2*(b[11]-pt0_Abd)^2)
  Vsp = 2 * vM_Splanchnic_vein * atan(consp*(b[11]-pt0_Abd)) / π

  conll = (π*C_Leg_vein)/(2*vM_Leg_vein)
  Cll = C_Leg_vein/(1+conll^2*(b[13]-pt0_Leg)^2)
  Vll = 2 * vM_Leg_vein * atan(conll*(b[13]-pt0_Leg)) / π

  conab = (π*C_Abd_veins)/(2*vM_Abd_vein)
  Cab = C_Abd_veins/(1+conab^2*(b[14]-pt0_Abd)^2)
  Vab = 2 * vM_Abd_vein * atan(conab*(b[14]-pt0_Abd)) / π
  A = zeros(N_STATE, N_STATE)
  F = zeros(N_STATE)

  F[1] = TBV -
    (v0_Asc_A + v0_BC_A + v0_UpBd_art + v0_Thor_A +
     v0_Abd_A + v0_Renal_art + v0_Splanchnic_art + v0_Leg_art +
     v0_UpBd_vein + v0_SVC + v0_Renal_vein + v0_Splanchnic_vein +
     v0_Leg_vein + v0_Abd_veins + v0_Thor_IVC + v0_ra +
     v0_rv + v0pa + v0pv + v0_la +
     v0_lv) +
    (1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+
     Cpv+C_Thor_IVC+C_SVC+C_Asc_A+C_BC_A+
     C_Thor_A)*pₜₕ+
     +pt0_UB*(C_UpBd_art+C_UpBd_vein) +
     pt0_Thor*(C_Thor_A+C_SVC+C_Asc_A+C_Thor_IVC+C_BC_A+1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+Cpv) +
     pt0_Abd*(C_Abd_A+C_Renal_art+C_Renal_vein+C_Splanchnic_art) +
     pt0_Leg*(C_Leg_art)-
    (C_UpBd_vein*b[4]+C_Renal_vein*b[9]+Vsp+Vll+Vab+C_Thor_IVC*b[15]+
     C_SVC*b[5]+1/Era0*b[16]+1/Ed_rv*b[17]+
     Cpa*b[19]+Cpv*b[20]+1/Ela0*b[21]+
     1/Ed_lv*b[22]+C_Asc_A*b[1]+C_BC_A*b[2]+
     C_UpBd_art*b[3]+C_Thor_A*b[6]+C_Abd_A*b[7]+
     C_Renal_art*b[8]+C_Splanchnic_art*b[10]+C_Leg_art*b[12])

  F[2] = (-1/R_BC_A-1/R_Thor_A-Tsys/T/R_Asc_A)*b[1] +
    b[2]/R_BC_A+b[6]/R_Thor_A+Tsys/T/R_Asc_A*b[23]

  F[3] = (-1/R_BC_A-1/R_Thor_A)*b[2] + b[1]/R_BC_A +
    b[3]/R_Thor_A

  F[4] = (-1/R_Thor_A-1/R_UpBd_cap)*b[3] + b[2]/R_Thor_A +
    b[4]/R_UpBd_cap

  F[5] = (-1/R_UpBd_cap-1/R_UpBd_vein)*b[4] + b[3]/R_UpBd_cap +
    b[4]/R_UpBd_vein

  F[6] = (-1/R_Thor_A-1/R_Abd_A)*b[6] + b[1]/R_Thor_A +
    b[7]/R_Abd_A

  F[7] = b[8]/R_Renal_art + b[10]/R_Splanchnic_art + b[12]/R_Leg_art -
    (1/R_Abd_A+1/R_Renal_art+1/R_Splanchnic_art+1/R_Leg_art)*
    b[7] + b[6]/R_Abd_A

  F[8] = (-1/R_Renal_art-1/R_Renal_cap)*b[8] + b[7]/R_Renal_art +
    b[9]/R_Renal_cap

  F[9] = (-1/R_Renal_cap-1/R_Renal_vein)*b[9] + b[8]/R_Renal_cap +
    b[14]/R_Renal_vein

  F[10] = (-1/R_Splanchnic_art-1/R_Splanchnic_cap)*b[10] + b[7]/R_Splanchnic_art +
    b[11]/R_Splanchnic_cap

  F[11] = (-1/R_Splanchnic_cap-1/R_Splanchnic_vein)*b[11] + b[10]/R_Splanchnic_cap +
    b[14]/R_Splanchnic_vein

  F[12] = (-1/R_Leg_art-1/R_Leg_cap)*b[12] + b[7]/R_Leg_art +
    b[13]/R_Leg_cap

  F[13] = (-1/R_Leg_cap-1/R_Leg_vein)*b[13] + b[12]/R_Leg_cap +
    b[14]/R_Leg_vein

  F[14] = b[13]/R_Leg_vein + b[11]/R_Splanchnic_vein + b[9]/R_Renal_vein -
    (1/R_Renal_vein+1/R_Splanchnic_vein+1/R_Leg_vein+1/R_Abd_veins)*b[14]
    + b[15]/R_Abd_veins

  F[15] = (-1/R_Abd_veins-1/R_Thor_IVC)*b[15] + b[14]/R_Abd_veins +
    b[16]/R_Thor_IVC

  F[16] = (-Tdias/T/R_tv-1/R_Thor_IVC-1/R_SVC)*b[16] +
    b[5]/R_SVC + b[15]/R_Thor_IVC + Tdias/T/R_tv*b[17]

  F[17] = -Tdias/T/R_tv*(b[17]-b[16]) - Tsys/T/Rpa *
    (b[18]-b[19])

  F[18] = (1/Ees_rv+Tsys/Rpa)*b[18]-1/Ees_rv*pₜₕ-
    1/Ed_rv*b[17]-Tsys/Rpa*b[19]+pₜₕ*1/Ed_rv

  F[19] = (1/Rpc+Tsys/T/Rpa)*b[19] - Tsys/T/Rpa*
    b[18] - b[20]/Rpc

  F[20] = (-1/Rpc-1/Rpv)*b[20] + b[19]/Rpc +
    b[21]/Rpv

  F[21] = (-1/Rpv-Tdias/T/R_mv)*b[21]+b[20]/Rpv+
    b[22]*Tdias/T/R_mv

  F[22] = -Tdias/T/R_mv*(b[21]-b[22])+Tsys/T/R_Asc_A*
    (b[23]-b[1])

  F[23] = (1/Ees_lv+Tsys/R_Asc_A)*b[23]-1/Ees_lv*pₜₕ-
    Tsys/R_Asc_A*b[1]-1/Ed_lv*(b[22]-pₜₕ)

    A[1,1]  = C_Asc_A          # Aortic arch
    A[1,2]  = C_BC_A          # Common Carotids
    A[1,3]  = C_UpBd_art         # Cerebral Arteries
    A[1,4]  = C_UpBd_vein           # Cerebral veins
    A[1,5]  = C_SVC           # Jugular veins
    A[1,6]  = C_Thor_A           # Thoracic aorta
    A[1,7]  = C_Abd_A          # Abdominal aorta
    A[1,8]  = C_Renal_art          # Renal arteries
    A[1,9]  = C_Renal_vein           # Renal veins
    A[1,10] = C_Splanchnic_art          # Splanchnic arteries
    A[1,11] = Csp           # Splanchnic veins
    A[1,12] = C_Leg_art          # Leg arteries
    A[1,13] = Cll           # Leg veins
    A[1,14] = Cab           # Abdominal veins
    A[1,15] = C_Thor_IVC           # Vena cava inf.
    A[1,16] = 1/Era0           # Right atrium
    A[1,17] = 1/Ed_rv          # Right ventricle diastole
    A[1,18] = 0.0                 # Right ventricle systole
    A[1,19] = Cpa           # Pulmonary arteries
    A[1,20] = Cpv           # Pulmonary veins
    A[1,21] = 1/Ela0           # Left atrium
    A[1,22] = 1/Ed_lv           # Left ventricle diastole
    A[1,23] = 0.0                 # Left ventricle systole

    # Left ventricular outflow = aortic arch outflow
    A[2,1]  = 1/R_BC_A + 1/R_Thor_A + Tsys/T/R_Asc_A
    A[2,2]  = -1/R_BC_A
    A[2,6]  = -1/R_Thor_A
    A[2,23] = -Tsys/T/R_Asc_A

    # Common carotid arterial flow = cerebral arterial flow
    A[3,1] = -1/R_BC_A
    A[3,2] = 1/R_BC_A + 1/R_UpBd_art
    A[3,3] = -1/R_UpBd_art

    # Upper body microflow.
    A[4,2] = -1/R_UpBd_art
    A[4,3] = 1/R_UpBd_art + 1/R_UpBd_cap
    A[4,4] = -1/R_UpBd_cap


    # Upper body veins to superior vena cava
    A[5,3]  = -1/R_UpBd_cap
    A[5,4]  = 1/R_UpBd_cap + 1/R_UpBd_vein
    A[5,5]  = -1/R_UpBd_vein

    # Thoracic aorta to abdominal aorta
    A[6,1] = -1/R_Thor_A
    A[6,6] = 1/R_Thor_A + 1/R_Abd_A
    A[6,7] = -1/R_Abd_A

    # Abdominal aorta to systemics arteries
    A[7,6]  = -1/R_Abd_A
    A[7,7]  = 1/R_Abd_A + 1/R_Renal_art + 1/R_Splanchnic_art + 1/R_Leg_art
    A[7,8]  = -1/R_Renal_art
    A[7,10] = -1/R_Splanchnic_art
    A[7,12] = -1/R_Leg_art

    # Abdominal aorta to renal vein
    A[8,7] = -1/R_Renal_art
    A[8,8] = 1/R_Renal_art + 1/R_Renal_cap
    A[8,9] = -1/R_Renal_cap

    # Renal artery to abdominal veins
    A[9,8]  = -1/R_Renal_cap
    A[9,9]  = 1/R_Renal_cap + 1/R_Renal_vein
    A[9,14] = -1/R_Renal_vein

    # Abdominal aorta to splanchnic veins
    A[10,7]  = -1/R_Splanchnic_art
    A[10,10] = 1/R_Splanchnic_art + 1/R_Splanchnic_cap
    A[10,11] = -1/R_Splanchnic_cap

    # Flow through splanchnic compartment.
    A[11,10] = -1/R_Splanchnic_cap
    A[11,11] = 1/R_Splanchnic_cap + 1/R_Splanchnic_vein
    A[11,14] = -1/R_Splanchnic_vein

    # Abdominal aorta to leg arteries
    A[12,7]  = -1/R_Leg_art
    A[12,12] = 1/R_Leg_art + 1/R_Leg_cap
    A[12,13] = -1/R_Leg_cap

    # Leg arteries to abdominal veins
    A[13,12] = -1/R_Leg_cap
    A[13,13] = 1/R_Leg_cap + 1/R_Leg_vein
    A[13,14] = -1/R_Leg_vein

    # Parallel venous outflow to abdominal veins
    A[14,9]  = -1/R_Renal_vein
    A[14,11] = -1/R_Splanchnic_vein
    A[14,13] = -1/R_Leg_vein
    A[14,14] = 1/R_Renal_vein + 1/R_Splanchnic_vein + 1/R_Leg_vein + 1/R_Abd_veins
    A[14,15] = -1/R_Abd_veins

    # Abdominal veins to right atrium
    A[15,14] = -1/R_Abd_veins
    A[15,15] = 1/R_Abd_veins + 1/R_Thor_IVC
    A[15,16] = -1/R_Thor_IVC

    # Right atrial inflow to right ventricle
    A[16,5]  = -1/R_SVC
    A[16,15] = -1/R_Thor_IVC
    A[16,16] = Tdias/T/R_tv + 1/R_Thor_IVC + 1/R_SVC
    A[16,17] = -Tdias/T/R_tv

    # Tricuspid flow to pulmonary artery
    A[17,16] = -Tdias/T/R_tv
    A[17,17] = Tdias/T/R_tv
    A[17,18] = Tsys/T/Rpa
    A[17,19] = -Tsys/T/Rpa

    # Left atrial inflow.
    A[18,17] = 1/Ed_rv
    A[18,18] = -1/Ees_rv - Tsys/Rpa
    A[18,19] = Tsys/Rpa

    # Left atrial inflow.
    A[19,18]  =Tsys/T/Rpa
    A[19,19]  =-1/Rpc-Tsys/T/Rpa
    A[19,20]  =1/Rpc

    # Left atrial inflow.
    A[20,19]  =-1/Rpc
    A[20,20]  = 1/Rpc+1/Rpv
    A[20,21]  =-1/Rpv

    # Left atrial inflow.
    A[21,20]  =-1/Rpv
    A[21,21]  = 1/Rpv+Tdias/T/R_mv
    A[21,22]  =-Tdias/T/R_mv

    # Left atrial inflow.
    A[22,1]   = Tsys/T/R_Asc_A
    A[22,21]  = Tdias/T/R_mv
    A[22,22]  =-Tdias/T/R_mv
    A[22,23]  =-Tsys/T/R_Asc_A

    # Left atrial inflow.
    A[23,1]   = Tsys/R_Asc_A
    A[23,22]  = 1/Ed_lv
    A[23,23]  =-1/Ees_lv-Tsys/R_Asc_A

   F1 = A \ F

    b .+= F1

    consp = (π*C_Splanchnic_vein)/(2*vM_Splanchnic_vein)
    Csp = C_Splanchnic_vein/(1+consp^2*(b[11]-pt0_Abd)^2)
    Vsp = 2 * vM_Splanchnic_vein * atan(consp*(b[11]-pt0_Abd)) / π

    conll = (π*C_Leg_vein)/(2*vM_Leg_vein)
    Cll = C_Leg_vein/(1+conll^2*(b[13]-pt0_Leg)^2)
    Vll = 2 * vM_Leg_vein * atan(conll*(b[13]-pt0_Leg)) / π

    conab = (π*C_Abd_veins)/(2*vM_Abd_vein)
    Cab = C_Abd_veins/(1+conab^2*(b[14]-pt0_Abd)^2)
    Vab = 2 * vM_Abd_vein * atan(conab*(b[14]-pt0_Abd)) / π

    vol_error = TBV -
    (v0_Asc_A + v0_BC_A + v0_UpBd_art + v0_Thor_A +
     v0_Abd_A + v0_Renal_art + v0_Splanchnic_art + v0_Leg_art +
     v0_UpBd_vein + v0_SVC + v0_Renal_vein + v0_Splanchnic_vein +
     v0_Leg_vein + v0_Abd_veins + v0_Thor_IVC + v0_ra +
     v0_rv + v0pa + v0pv + v0_la +
     v0_lv) +
    (1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+
     Cpv+C_Thor_IVC+C_SVC+C_Asc_A+C_BC_A+
     C_Thor_A)*pₜₕ+
     +pt0_UB*(C_UpBd_art+C_UpBd_vein) +
     pt0_Thor*(C_Thor_A+C_SVC+C_Asc_A+C_Thor_IVC+C_BC_A+1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+Cpv) +
     pt0_Abd*(C_Abd_A+C_Renal_art+C_Renal_vein+C_Splanchnic_art) +
     pt0_Leg*(C_Leg_art)-
    (C_UpBd_vein*b[4]+C_Renal_vein*b[9]+Vsp+Vll+Vab+C_Thor_IVC*b[15]+
     C_SVC*b[5]+1/Era0*b[16]+1/Ed_rv*b[17]+
     Cpa*b[19]+Cpv*b[20]+1/Ela0*b[21]+
     1/Ed_lv*b[22]+C_Asc_A*b[1]+C_BC_A*b[2]+
     C_UpBd_art*b[3]+C_Thor_A*b[6]+C_Abd_A*b[7]+
     C_Renal_art*b[8]+C_Splanchnic_art*b[10]+C_Leg_art*b[12])

  return vol_error
end

"""
Create the IC matrices
This code section defines the A and b matrices for the initial guess at the compartment pressures using the zero pressure filling volumes of the nonlinear compartments. The mnewt!() function is then called to iteratively solve for the true pressures. Routine originally defined by Heldt (2004) and modified by Whittle (2023) to include the external tissue pressures.
"""

# Linearized volume equation.
A[1,1]  = C_Asc_A          # Aortic arch
A[1,2]  = C_BC_A          # Common Carotids
A[1,3]  = C_UpBd_art         # Cerebral Arteries
A[1,4]  = C_UpBd_vein           # Cerebral veins
A[1,5]  = C_SVC           # Jugular veins
A[1,6]  = C_Thor_A           # Thoracic aorta
A[1,7]  = C_Abd_A          # Abdominal aorta
A[1,8]  = C_Renal_art          # Renal arteries
A[1,9]  = C_Renal_vein           # Renal veins
A[1,10] = C_Splanchnic_art          # Splanchnic arteries
A[1,11] = C_Splanchnic_vein           # Splanchnic veins
A[1,12] = C_Leg_art          # Leg arteries
A[1,13] = C_Leg_vein           # Leg veins
A[1,14] = C_Abd_veins           # Abdominal veins
A[1,15] = C_Thor_IVC           # Vena cava inf.
A[1,16] = 1/Era0           # Right atrium
A[1,17] = 1/Ed_rv          # Right ventricle diastole
A[1,18] = 0.0                 # Right ventricle systole
A[1,19] = Cpa           # Pulmonary arteries
A[1,20] = Cpv           # Pulmonary veins
A[1,21] = 1/Ela0           # Left atrium
A[1,22] = 1/Ed_lv           # Left ventricle diastole
A[1,23] = 0.0                 # Left ventricle systole

# Left ventricular outflow = aortic arch outflow
A[2,1]  = 1/R_BC_A + 1/R_Thor_A + Tsys/T/R_Asc_A
A[2,2]  = -1/R_BC_A
A[2,6]  = -1/R_Thor_A
A[2,23] = -Tsys/T/R_Asc_A

# Common carotid arterial flow = cerebral arterial flow
A[3,1] = -1/R_BC_A
A[3,2] = 1/R_BC_A + 1/R_UpBd_art
A[3,3] = -1/R_UpBd_art

# Upper body microflow.
A[4,2] = -1/R_UpBd_art
A[4,3] = 1/R_UpBd_art + 1/R_UpBd_cap
A[4,4] = -1/R_UpBd_cap


# Upper body veins to superior vena cava
A[5,3]  = -1/R_UpBd_cap
A[5,4]  = 1/R_UpBd_cap + 1/R_UpBd_vein
A[5,5]  = -1/R_UpBd_vein

# Thoracic aorta to abdominal aorta
A[6,1] = -1/R_Thor_A
A[6,6] = 1/R_Thor_A + 1/R_Abd_A
A[6,7] = -1/R_Abd_A

# Abdominal aorta to systemics arteries
A[7,6]  = -1/R_Abd_A
A[7,7]  = 1/R_Abd_A + 1/R_Renal_art + 1/R_Splanchnic_art + 1/R_Leg_art
A[7,8]  = -1/R_Renal_art
A[7,10] = -1/R_Splanchnic_art
A[7,12] = -1/R_Leg_art

# Abdominal aorta to renal vein
A[8,7] = -1/R_Renal_art
A[8,8] = 1/R_Renal_art + 1/R_Renal_cap
A[8,9] = -1/R_Renal_cap

# Renal artery to abdominal veins
A[9,8]  = -1/R_Renal_cap
A[9,9]  = 1/R_Renal_cap + 1/R_Renal_vein
A[9,14] = -1/R_Renal_vein

# Abdominal aorta to splanchnic veins
A[10,7]  = -1/R_Splanchnic_art
A[10,10] = 1/R_Splanchnic_art + 1/R_Splanchnic_cap
A[10,11] = -1/R_Splanchnic_cap

# Flow through splanchnic compartment.
A[11,10] = -1/R_Splanchnic_cap
A[11,11] = 1/R_Splanchnic_cap + 1/R_Splanchnic_vein
A[11,14] = -1/R_Splanchnic_vein

# Abdominal aorta to leg arteries
A[12,7]  = -1/R_Leg_art
A[12,12] = 1/R_Leg_art + 1/R_Leg_cap
A[12,13] = -1/R_Leg_cap

# Leg arteries to abdominal veins
A[13,12] = -1/R_Leg_cap
A[13,13] = 1/R_Leg_cap + 1/R_Leg_vein
A[13,14] = -1/R_Leg_vein

# Parallel venous outflow to abdominal veins
A[14,9]  = -1/R_Renal_vein
A[14,11] = -1/R_Splanchnic_vein
A[14,13] = -1/R_Leg_vein
A[14,14] = 1/R_Renal_vein + 1/R_Splanchnic_vein + 1/R_Leg_vein + 1/R_Abd_veins
A[14,15] = -1/R_Abd_veins

# Abdominal veins to right atrium
A[15,14] = -1/R_Abd_veins
A[15,15] = 1/R_Abd_veins + 1/R_Thor_IVC
A[15,16] = -1/R_Thor_IVC

# Right atrial inflow to right ventricle
A[16,5]  = -1/R_SVC
A[16,15] = -1/R_Thor_IVC
A[16,16] = Tdias/T/R_tv + 1/R_Thor_IVC + 1/R_SVC
A[16,17] = -Tdias/T/R_tv

# Tricuspid flow to pulmonary artery
A[17,16] = -Tdias/T/R_tv
A[17,17] = Tdias/T/R_tv
A[17,18] = Tsys/T/Rpa
A[17,19] = -Tsys/T/Rpa

# Left atrial inflow.
A[18,17] = 1/Ed_rv
A[18,18] = -1/Ees_rv - Tsys/Rpa
A[18,19] = Tsys/Rpa

# Left atrial inflow.
A[19,18]  =-Tsys/T/Rpa
A[19,19]  = 1/Rpc+Tsys/T/Rpa
A[19,20]  =-1/Rpc

# Left atrial inflow.
A[20,19]  =-1/Rpc
A[20,20]  = 1/Rpc+1/Rpv
A[20,21]  =-1/Rpv

# Left atrial inflow.
A[21,20]  =-1/Rpv
A[21,21]  = 1/Rpv+Tdias/T/R_mv
A[21,22]  =-Tdias/T/R_mv

# Left atrial inflow.
A[22,1]   = Tsys/T/R_Asc_A
A[22,21]  = Tdias/T/R_mv
A[22,22]  =-Tdias/T/R_mv
A[22,23]  =-Tsys/T/R_Asc_A

# Left atrial inflow.
A[23,1]   = Tsys/R_Asc_A
A[23,22]  = 1/Ed_lv
A[23,23]  =-1/Ees_lv-Tsys/R_Asc_A

b[1]  = TBV -
  (v0_Asc_A+v0_BC_A+v0_UpBd_art+v0_Thor_A +
    v0_Abd_A+v0_Renal_art+v0_Splanchnic_art+v0_Leg_art +
    v0_UpBd_vein+v0_SVC+v0_Renal_vein+v0_Splanchnic_vein +
    v0_Leg_vein+v0_Abd_veins+v0_Thor_IVC +
    v0_ra+v0_rv+v0pa+v0pv+v0_la+v0_lv) +
    pₜₕ*(1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+
      Cpa+Cpv+C_Thor_IVC+C_SVC+
      C_Asc_A+C_BC_A+C_Thor_A) +
    pt0_UB*(C_UpBd_art+C_UpBd_vein) +
    pt0_Thor*(C_Thor_A+C_SVC+C_Asc_A+C_Thor_IVC+C_BC_A+1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+Cpv) +
    pt0_Abd*(C_Abd_A+C_Renal_art+C_Renal_vein+C_Splanchnic_art+
      C_Splanchnic_vein+C_Abd_veins) +
    pt0_Leg*(C_Leg_art+C_Leg_vein)

b[18] = pₜₕ*(1/Ed_rv-1/Ees_rv)
b[23] = pₜₕ*(1/Ed_lv-1/Ees_lv)

"""
Solve to find the initial pressure vector
In the Heldt C implementation, this was a long routine. Here in Julia, the initial guess is a single line. A short function solve_volume!() is defined to iterate using the previously defined mnewt!() function. This converges on the true inital pressure vector.
"""

x = A \ b
volume_error = abs(mnewt!(x))
function solve_volume!(b; tol = 1e-11)
  volume_error = abs(mnewt!(b))
  while volume_error > tol
      volume_error = abs(mnewt!(b))
  end
  return b
end

x = solve_volume!(x)

"""
Debug Script
This short code section is used as a sanity check to ensure that the total blood volume calculated from the initial compartment pressures is equal to the total blood volume defined in the model parameters file.
"""

consp = (π*C_Splanchnic_vein)/(2*vM_Splanchnic_vein)
Csp = C_Splanchnic_vein/(1+consp^2*(x[11]-pt0_Abd)^2)
Vsp = 2 * vM_Splanchnic_vein * atan(consp*(x[11]-pt0_Abd)) / π

conll = (π*C_Leg_vein)/(2*vM_Leg_vein)
Cll = C_Leg_vein/(1+conll^2*(x[13]-pt0_Leg)^2)
Vll = 2 * vM_Leg_vein * atan(conll*(x[13]-pt0_Leg)) / π

conab = (π*C_Abd_veins)/(2*vM_Abd_vein)
Cab = C_Abd_veins/(1+conab^2*(x[14]-pt0_Abd)^2)
Vab = 2 * vM_Abd_vein * atan(conab*(x[14]-pt0_Abd)) / π

TBV -
    (v0_Asc_A + v0_BC_A + v0_UpBd_art + v0_Thor_A +
     v0_Abd_A + v0_Renal_art + v0_Splanchnic_art + v0_Leg_art +
     v0_UpBd_vein + v0_SVC + v0_Renal_vein + v0_Splanchnic_vein +
     v0_Leg_vein + v0_Abd_veins + v0_Thor_IVC + v0_ra +
     v0_rv + v0pa + v0pv + v0_la +
     v0_lv) +
    (1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+
     Cpv+C_Thor_IVC+C_SVC+C_Asc_A+C_BC_A+
     C_Thor_A)*pₜₕ+
     +pt0_UB*(C_UpBd_art+C_UpBd_vein) +
     pt0_Thor*(C_Thor_A+C_SVC+C_Asc_A+C_Thor_IVC+C_BC_A+1/Ela0+1/Era0+1/Ed_lv+1/Ed_rv+Cpa+Cpv) +
     pt0_Abd*(C_Abd_A+C_Renal_art+C_Renal_vein+C_Splanchnic_art) +
     pt0_Leg*(C_Leg_art)-
    (C_UpBd_vein*x[4]+C_Renal_vein*x[9]+Vsp+Vll+Vab+C_Thor_IVC*x[15]+
     C_SVC*x[5]+1/Era0*x[16]+1/Ed_rv*x[17]+
     Cpa*x[19]+Cpv*x[20]+1/Ela0*x[21]+
     1/Ed_lv*x[22]+C_Asc_A*x[1]+C_BC_A*x[2]+
     C_UpBd_art*x[3]+C_Thor_A*x[6]+C_Abd_A*x[7]+
     C_Renal_art*x[8]+C_Splanchnic_art*x[10]+C_Leg_art*x[12])