"""
This file contains the model parameters for the cardiovascular system simulation. The parameters are organized into different sections, including system parameters, segment pressures, heart parameters, pulmonary system parameters, microvascular resistances, arterial system parameters, venous system parameters, interstitial compartment parameters, carotid sinus parameters, and reflex system parameters.
"""

module ModelParams

"""
Global System Parameters
These parameters define global anthropometric and physiological constants for the model.
"""
HRₙₒₘ = 67.0  # Nominal Heart Rate (bpm)
TBV = 5300.0 # Total Blood Volume (ml)
ρ_b = 1060.0 # Blood Density (kg/m^3)
Pa2mmHg = 0.00750062 # Conversion factor from Pascal to mmHg
mmHg2dynecm2 = 10/Pa2mmHg # Conversion factor from mmHg to dyne/cm^2
cmH2O2mmHg = 0.73555912 # Conversion factor from cmH2O to mmHg
η_b = 0.004 # Blood Viscosity (Pa.s)

p₀ = 0.0 # Atmospheric pressure (mmHg)
pₚₗ₀ = -5.0 * cmH2O2mmHg # Initial Pleural pressure (mmHg)

p_col = 5.0 # Vascular Collapsing Pressure (negative) (mmHg)

"""
Segment Pressures
These parameters define the external static pressures on different segments of the cardiovascular system.
"""
p_abd = 0.0 # Abdominal Pressure (mmHg)
p_icp = 10.0 # Intracranial Pressure (mmHg)

con_default = 2.0 # Default hydrostatic conversion factor

"""
Tissue Parameters
These parameters define the tissue properties, including the density of fat-free tissue and the radii of different segments of the body. The radii are used set the hydrostatic pressure from external tissue weight.
"""
ρ_fft = 1100.0  # Fat Free Tissue Density (kg/m^3)

rad_Head = 10.0 # Radius of the Head Segment (cm)
rad_Neck = 7.0 # Radius of the Neck Segment (cm)
rad_UB = 7.0 # Radius of the Brachiocephalic Segment (cm)
rad_Thor = 14.0 # Radius of the Thoracic Segment (cm)
rad_Abd = 12.0 # Radius of the Abdominal Segment (cm)
rad_Leg = 10.0 # Radius of the Leg Segment (cm)

"""
Heart Parameters
These parameters define the heart's properties, including systolic and elastance and zero pressure volumes for different chambers, as well as the timings of the atrial-ventricular cycle. The mitral and tricuspid valves are also defined here (aortic and pulmonary valves are lumped with their respective arteries).
"""

#### Cardiac Timings
RRₙₒₘ = 60.0 / HRₙₒₘ # Nominal RR interval (s)
τₐᵥ = 0.12 # PR interval (A-V delay) (s)
τₐₛ = 0.2 # Atrium systolic contraction time (s)
τᵥₛ = 0.3 # QT interval (Ventricle systolic contraction time) (s)

#### Right Atrium
v0_ra = 14.0 # Right Atrium zero pressure volume (ml)
Ed_ra = 0.30 # Right Atrium diastolic elastance (mmHg/ml)
Ees_ra = 0.74 # Right Atrium end systolic elastance (mmHg/ml)

#### Right Ventricle
v0_rv = 46.0 # Right Ventricle zero pressure volume (ml)
Ed_rv = 0.05 # Right Ventricle diastolic elastance (mmHg/ml)
Ees_rv = 1.3 # Right Ventricle end systolic elastance (mmHg/ml)
Elimit_rv = 3.0 # Right Ventricle elastance limit (mmHg/ml)

#### Left Atrium
v0_la = 24.0 # Left Atrium zero pressure volume (ml)
Ed_la = 0.50 # Left Atrium diastolic elastance (mmHg/ml)
Ees_la = 0.61 # Left Atrium end systolic elastance (mmHg/ml)

#### Left Ventricle
v0_lv = 55.0 # Left Ventricle zero pressure volume (ml)
Ed_lv = 0.11 # Left Ventricle diastolic elastance (mmHg/ml)
Ees_lv = 2.5 # Left Ventricle end systolic elastance (mmHg/ml)
Elimit_lv = 8.0 # Left Ventricle elastance limit (mmHg/ml)

#### Tricuspid Valve
Ann_tv = 8.0 # Annulus area of the tricuspid valve (cm^2)
Kvc_tv = 0.04 # Tricuspid valve closing rate coefficient (cm^2/(dynes*s))
Kvo_tv = 0.03 # Tricuspid valve opening rate coefficient (cm^2/(dynes*s))

#### Pulmonary Valve
Leff_pv = 1.6 # Effective length of the pulmonary valve (cm) Pulmonary trunk length
Ann_pv = 7.1 # Annulus area of the pulmonary valve (cm^2)
Kvc_pv = 0.02 # Pulmonary valve closing rate coefficient (cm^2/(dynes*s))
Kvo_pv = 0.02 # Pulmonary valve opening rate coefficient (cm^2/(dynes*s))

#### Mitral Valve
Ann_mv = 7.7 # Annulus area of the mitral valve (cm^2)
Kvc_mv = 0.04 # Mitral valve closing rate coefficient (cm^2/(dynes*s))
Kvo_mv = 0.03 # Mitral valve opening rate coefficient (cm^2/(dynes*s))

#### Aortic Valve
Leff_av = 1.8 # Effective length of the aortic valve (cm) Aortic arch length
Ann_av = 6.8 # Annulus area of the aortic valve (cm^2)
Kvc_av = 0.012 # Aortic valve closing rate coefficient (cm^2/(dynes*s))
Kvo_av = 0.0120 # Aortic valve opening rate coefficient (cm^2/(dynes*s))

"""
Pulmonary System Parameters
These parameters define the properties of the pulmonary system, including the pulmonary artery and vein resistances, compliances, and zero pressure volumes. The pulmonary capillary microvascular resistance is also defined here.
"""

#### Pulmonary Artery
Rpa = 0.006 # Pulmonary Artery resistance (PRU)
Cpa = 3.4 # Pulmonary Artery compliance (ml/mmHg)
v0pa = 160.0 # Pulmonary Artery zero pressure volume (ml)
# Lpa is included in the valve model

#### Pulmonary Capillaries
Rpc = 0.07 # Pulmonary Capillary microvascular resistance (PRU)
h_Lungs = 23.6 # Height of the lungs (cm)

#### Pulmonary Vein parameters
Rpv = 0.006 # Pulmonary Vein resistance (PRU)
Cpv = 9.0 # Pulmonary Vein compliance (ml/mmHg)
v0pv = 430.0 # Pulmonary Vein zero pressure volume (ml)

"""
Coronary System Parameters
These parameters define the properties of the coronary system, including the coronary artery and vein resistances, compliances, and zero pressure volumes. The coronary microvascular resistance is also defined here. Values were taken from Albanese (2016) and subracted from splanchnic values in the original Heldt model.
"""

#### Coronary Artery
Rca = 0.34 # Coronary Artery resistance (PRU)
Cca = 0.15 # Coronary Artery compliance (ml/mmHg)
v0ca = 24.0 # Coronary Artery zero pressure volume (ml)
Lca = 0.047 # Coronary Artery inductance (mmHg.s^2/ml)

#### Coronary Capillaries
Rcc = 19.71 # Coronary Capillary microvascular resistance (PRU)

#### Coronary Vein parameters
Rcv = 0.224 # Coronary Vein resistance (PRU)
Ccv = 2.5 # Coronary Vein compliance (ml/mmHg)
v0cv = 98.2 # Coronary Vein zero pressure volume (ml)

"""
Microvascular Resistances
These parameters define the baseline microvascular resistances for different compartments of the body, including the upper body, renal, splanchnic, and leg compartments.
"""

R_Head_cap = 8.1 # Head microvascular resistance (PRU)
R_UpBd_cap = 9.7 # Upper Body microvascular resistance (PRU)
R_Renal_cap = 4.7 # Renal microvascular resistance (PRU)
R_Splanchnic_cap = 3.1 # Splanchnic microvascular resistance (PRU) (Increased from 2.8 for coronary arteries)
R_Leg_cap = 4.0 # Leg microvascular resistance (PRU)

"""
Arterial System Parameters
These parameters define the properties of the arterial system, including the resistances, compliances, zero pressure volumes, and the hydrostatic vertical lengths of different arterial compartments.
"""

#### Compartment 1: Ascending Aorta
R_Asc_A = 0.007 # Ascending Aorta resistance (PRU)
C_Asc_A = 0.28 # Ascending Aorta compliance (ml/mmHg)
v0_Asc_A = 21.0 # Ascending Aorta zero pressure volume (ml)
h_Asc_A = -10.0 # Ascending Aorta length (cm)
# L_Asc_A is included in the valve model

#### Compartment 2: Bracehocephalic Arteries
R_BC_A = 0.003 # Bracehocephalic Artery resistance (PRU)
C_BC_A = 0.13 # Bracehocephalic Artery compliance (ml/mmHg)
v0_BC_A = 5.0 # Bracehocephalic Artery zero pressure volume (ml)
h_BC_A = -4.5 # Bracehocephalic Artery length (cm)
L_BC_A = 0.00008 # Bracehocephalic Artery inductance (mmHg.s^2/ml)

#### Compartment 3: Upper Body Arteries
R_UpBd_art = 0.014 # Upper Body Artery resistance (PRU)
C_UpBd_art = 0.26 # Upper Body Artery compliance (ml/mmHg)
v0_UpBd_art = 72.0 # Upper Body Artery zero pressure volume (ml)
h_UpBd_art = 66.0 # Upper Body Artery length (cm)
con_UpBd_art = 3.0 # Upper Body Artery hydrostatic conversion
L_UpBd_art = 0.00141 # Upper Body Artery inductance (mmHg.s^2/ml)

#### Compartment 6: Thoracic Aorta
R_Thor_A = 0.011 # Thoracic Aorta resistance (PRU)
C_Thor_A = 0.21 # Thoracic Aorta compliance (ml/mmHg)
v0_Thor_A = 16.0 # Thoracic Aorta zero pressure volume (ml)
h_Thor_A = 16.0 # Thoracic Aorta length (cm)
L_Thor_A = 0.00016 # Thoracic Aorta inductance (mmHg.s^2/ml)

#### Compartment 7: Abdominal Aorta
R_Abd_A = 0.01 # Abdominal Aorta resistance (PRU)
C_Abd_A = 0.1 # Abdominal Aorta compliance (ml/mmHg)
v0_Abd_A = 10.0 # Abdominal Aorta zero pressure volume (ml)
h_Abd_A = 14.5 # Abdominal Aorta length (cm)
L_Abd_A = 0.00014 # Abdominal Aorta inductance (mmHg.s^2/ml)

#### Compartment 8: Renal Artery
R_Renal_art = 0.1 # Renal Artery resistance (PRU)
C_Renal_art = 0.21 # Renal Artery compliance (ml/mmHg)
v0_Renal_art = 20.0 # Renal Artery zero pressure volume (ml)
h_Renal_art = 0.0 # Renal Artery length (cm)
L_Renal_art = 0.00070 # Renal Artery inductance (mmHg.s^2/ml)

#### Compartment 10: Splanchnic Artery
R_Splanchnic_art = 0.07 # Splanchnic Artery resistance (PRU)
C_Splanchnic_art = 0.42 # Splanchnic Artery compliance (ml/mmHg)
v0_Splanchnic_art = 276.0 # Splanchnic Artery zero pressure volume (ml) (Subtracted coronary arteries)
h_Splanchnic_art = 10.0 # Splanchnic Artery length (cm)
L_Splanchnic_art = 0.00031 # Splanchnic Artery inductance (mmHg.s^2/ml)

#### Compartment 12: Leg Arteries
R_Leg_art = 0.09 # Leg Artery resistance (PRU)
C_Leg_art = 0.42 # Leg Artery compliance (ml/mmHg)
v0_Leg_art = 200.0 # Leg Artery zero pressure volume (ml)
h_Leg_art = 105.0 # Leg Artery length (cm)
con_Leg_art = 3.0 # Leg Artery hydrostatic conversion
L_Leg_art = 0.00277 # Leg Artery inductance (mmHg.s^2/ml)

#### Compartment H1: Common Carotid Arteries
R_CCA = 0.014 # Common Carotid Artery resistance (PRU)
C_CCA = 0.07 # Common Carotid Artery compliance (ml/mmHg)
v0_CCA = 20.0 # Common Carotid Artery zero pressure volume (ml)
h_CCA = -20.0 # Common Carotid Artery length (cm)
L_CCA = 0.0007 # Common Carotid Artery inductance (mmHg.s^2/ml)

#### Compartment H2: Head Arteries
R_Head_art = 0.014 # Head Artery resistance (PRU)
C_Head_art = 0.08 # Head Artery compliance (ml/mmHg)
v0_Head_art = 108.0 # Head Artery zero pressure volume (ml)
h_Head_art = -20.0 # Head Artery length (cm)
L_Head_art = 0.0007 # Head Artery inductance (mmHg.s^2/ml)

"""
Venous System Parameters
These parameters define the properties of the venous system, including the resistances, compliances, zero pressure volumes, and the hydrostatic vertical lengths of different venous compartments. For the nonlinear compartments, the maximum distending volumes are also defined.
"""

#### Compartment 4: Upper Body Veins
R_UpBd_vein = 0.11 # Upper Body Vein resistance (PRU)
C_UpBd_vein = 1.2 # Upper Body Vein compliance (ml/mmHg)
v0_UpBd_vein = 360.0 # Upper Body Veins zero pressure volume (ml)
h_UpBd_vein = -66.0 # Upper Body Veins length (cm)
con_UpBd_vein = 3.0 # Upper Body Veins hydrostatic conversion

#### Compartment 5: Superior Vena Cava
R_SVC = 0.028 # Superior Vena Cava resistance (PRU)
C_SVC = 1.3 # Superior Vena Cava compliance (ml/mmHg)
v0_SVC = 16.0 # Superior Vena Cava zero pressure volume (ml)
h_SVC = 14.5 # Superior Vena Cava length (cm)

#### Compartment 9: Renal Vein
R_Renal_vein = 0.11 # Renal Vein resistance (PRU)
C_Renal_vein = 5.0 # Renal Vein compliance (ml/mmHg)
v0_Renal_vein = 30.0 # Renal Vein zero pressure volume (ml)
h_Renal_vein = 0.0 # Renal Vein length (cm)
vₘᵢₙ_Renal_vein = 5.0 # The Renal vein has a minimum zero-pressure volume (Zamanian, 2007)

#### Compartment 11: Splanchnic Vein (nonlinear)
R_Splanchnic_vein = 0.07 # Splanchnic Vein resistance (PRU)
C_Splanchnic_vein = 50.0 # Splanchnic Vein compliance (ml/mmHg)
v0_Splanchnic_vein = 1047.8 # Splanchnic Vein zero pressure volume (ml) (Subtracted coronary veins)
h_Splanchnic_vein = -10.0 # Splanchnic Vein length (cm)
vM_Splanchnic_vein = 1565.0 # Splanchnic Vein maximum volume (ml)


#### Compartment 13: Leg Vein (nonlinear)
R_Leg_vein = 0.1 # Leg Vein resistance (PRU)
C_Leg_vein = 27.0 # Leg Vein compliance (ml/mmHg)
v0_Leg_vein = 716.0 # Leg Vein zero pressure volume (ml)
h_Leg_vein = -105.0 # Leg Vein length (cm)
vM_Leg_vein = 1043.0 # Leg Vein maximum volume (ml)
con_Leg_vein = 3.0 # Leg Vein hydrostatic conversion

#### Compartment 14: Abdominal Veins (nonlinear)
R_Abd_veins = 0.019 # Abdominal Veins resistance (PRU)
C_Abd_veins = 1.3 # Abdominal Veins compliance (ml/mmHg)
v0_Abd_veins = 79.0 # Abdominal Veins zero pressure volume (ml)
h_Abd_veins = -14.5 # Abdominal Veins length (cm)
vM_Abd_vein = 678.0 # Abdominal Veins maximum volume (ml)

#### Compartment 15: Thoracic IVC
R_Thor_IVC = 0.008 # Thoracic IVC resistance (PRU)
C_Thor_IVC = 0.5 # Thoracic IVC compliance (ml/mmHg)
v0_Thor_IVC = 33.0 # Thoracic IVC zero pressure volume (ml)
h_Thor_IVC = -6.0 # Thoracic IVC length (cm)

#### Compartment H3: Head Veins
R_Head_veins = 0.05 # Head Vein resistance (PRU)
C_Head_veins = 3.35 # Head Vein compliance (ml/mmHg)
v0_Head_veins = 250.0 # Head Vein zero pressure volume (ml)
h_Head_veins = 20.0 # Head Vein length (cm)

#### Compartment H4: Jugular Veins
R_Jugular_vein = 0.05 # Jugular Vein resistance (PRU)
C_Jugular_vein = 2.45 # Jugular Vein compliance (ml/mmHg)
v0_Jugular_vein = 35.0 # Jugular Vein zero pressure volume (ml)
h_Jugular_vein = 20.0 # Jugular Vein length (cm)

#### Vertebral Plexus
Rᵥₚ = 0.068 # Vertebral Plexus resistance (PRU)

"""
Interstitial Compartment Parameters
These parameters define the properties of the interstitial compartment, including the time constant and the maximum instantaneous interstitial volume in both tilt and LBNP. We also include the interstitial flow rates for the splanchnic, leg, and abdominal veins.
"""

τint = 276 # Interstitial time constant (s)
Vmax_tilt = 700 # Interstitial maximum volume in tilt (ml)
Vmax_lbnp = 491 # Interstitial maximum volume in LBNP (ml)

Phtotal_Qint = (abs(h_Splanchnic_vein) / con_default) + (abs(h_Leg_vein) / con_Leg_vein) + (abs(h_Abd_veins) / con_default)
Flow_Splanchnic_vein = (abs(h_Splanchnic_vein) / con_default) / Phtotal_Qint
Flow_Leg_vein = (abs(h_Leg_vein) / con_Leg_vein) / Phtotal_Qint
Flow_Abd_veins = (abs(h_Abd_veins) / con_default) / Phtotal_Qint

"""
Carotid Sinus
The carotid sinus offset above the brachiocephalic arteries is defined here. This is used to calculate the hydrostatic pressure difference between the carotid sinus and the heart for the arterial baroreflex.
"""

h_cs = -14.0 # Carotid Sinus Offset (cm)

"""
Heldt Reflex System Parameters
These parameters define the properties of the arterial baroreflex (ABR) and cardiopulmonary reflex (CPR). The afferent arms include a state space filter with  matrices A, B, and C, as well as the set points and afferent gains. The transfer functions are defined by a delay, a peak, and an end timing and are dynamically set to have unit area impulse response. The delay is defined by a Pade approximation, the order of the approximation can be set here to give a sharper response (more computationally intensive). Finally, the static gains for the ABR and CPR are also defined here.
"""

######################
# Afferent Arms
######################

#### State Space Filter

A_9x9 = [
  -9.45889157397485 -116.684730988048 -716.776540007884 -4286.18041417701 -16964.9739125645 -56703.9505500527 -128365.488100833 -200510.747016871 -151338.279894609
  1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
] # A matrix for the state space filter
B_9x1 = [1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0][:, :] # B matrix for the state space filter
C_1x9 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 151338.279894609] # C matrix for the state space filter
x0 = ones(9) # Initial condition for the state space filter

#### Afferent Set Points and Gains

p_cpr = 6.0 # CPR Set Point (mmHg)
gain_cpr = 5.0 # CPR Afferent Gain

p_abr = 91.0 # ABR Set Point (mmHg)
gain_abr = 18.0 # ABR Afferent Gain

######################
# Transfer Functions
######################

#### Pade Approximation

reflex_delay_order = 2 # Order of the Pade delay
reflex_delay_init = zeros(reflex_delay_order) # Initial condition for the delay

#### Transfer Function Timings

abr_αr_delay = 2.5 # Arterial Baroreflex α-sympathetic resistance delay (s)
abr_αr_peak = 3.5 # Arterial Baroreflex α-sympathetic resistance peak (s)
abr_αr_end = 30.0 # Arterial Baroreflex α-sympathetic resistance end (s)

abr_αv_delay = 5.0 # Arterial Baroreflex α-sympathetic volume delay (s)
abr_αv_peak = 10.0 # Arterial Baroreflex α-sympathetic volume peak (s)
abr_αv_end = 42.0 # Arterial Baroreflex α-sympathetic volume end (s)

abr_β_delay = 2.5 # Arterial Baroreflex β-sympathetic delay (s)
abr_β_peak = 3.5 # Arterial Baroreflex β-sympathetic peak (s)
abr_β_end = 15.0 # Arterial Baroreflex β-sympathetic end (s)

abr_para_delay = 0.59 # Arterial Baroreflex parasympathetic delay (s)
abr_para_peak = 0.70 # Arterial Baroreflex parasympathetic peak (s)
abr_para_end = 1.0 # Arterial Baroreflex parasympathetic end (s)

cpr_αr_delay = 2.5 # Cardiopulmonary Reflex α-sympathetic resistance delay (s)
cpr_αr_peak = 5.5 # Cardiopulmonary Reflex α-sympathetic resistance peak (s)
cpr_αr_end = 35.0 # Cardiopulmonary Reflex α-sympathetic resistance end (s)

cpr_αv_delay = 5.0 # Cardiopulmonary Reflex α-sympathetic volume delay (s)
cpr_αv_peak = 9.0 # Cardiopulmonary Reflex α-sympathetic volume peak (s)
cpr_αv_end = 40.0 # Cardiopulmonary Reflex α-sympathetic volume end (s)

######################
# Efferent Static Gains
######################

#### Arterial Baroreflex

Gabr_r = -0.05 # ABR Resistance Gain (PRU/mmHg) (Previous -0.13)

Gabr_vub = 5.0 # ABR Upper Body Volume Gain (ml/mmHg) (Previous 5.3)
Gabr_vre = 2.0 # ABR Renal Volume Gain (ml/mmHg) (Previous 1.3)
Gabr_vsp = 13.0 # ABR Splanchnic Volume Gain (ml/mmHg) (Previous 13.3)
Gabr_vlb = 7.0 # ABR Lower Body Volume Gain (ml/mmHg) (Previous 6.7)

Gabr_erv = -0.037 # ABR RV Elastance Gain (1/ml) (Equivalent to contractility gain of 0.022 ml/mmHg^2) (Previous Contractility Gain 0.021)
Gabr_elv = -0.044 # ABR LV Elastance Gain (1/ml) (Equivalent to contractility gain of 0.007 ml/mmHg^2) (Previous Contractility Gain 0.014)

Gabr_rrsymp = 0.009 # ABR RR Interval Sympathetic Gain (s/mmHg) (Previous 0.012)
Gabr_rrpara = 0.009 # ABR RR Interval Parasympathetic Gain (s/mmHg)

#### Cardiopulmonary Reflex

Gcpr_r = -0.05 # CPR Resistance Gain (PRU/mmHg) (Previous -0.3)

Gcpr_vub = 13.0 # CPR Upper Body Volume Gain (ml/mmHg) (Previous 13.5)
Gcpr_vre = 3.0 # CPR Renal Volume Gain (ml/mmHg) (Previous 2.7)
Gcpr_vsp = 64.0 # CPR Splanchnic Volume Gain (ml/mmHg)
Gcpr_vlb = 30.0 # CPR Lower Body Volume Gain (ml/mmHg)

"""
Lung Model
"""

#### Respiratory Muscle Parameters

p_musmin = -5.0 * cmH2O2mmHg # Minimum respiratory muscle pressure (mmHg)
RespRateₙₒₘ = 12.0 # Breathing Rate (breaths/min)
IE_ratio = 0.6 # Inspiratory to Expiratory Ratio

#### Lung Parameters

Rml = 1.021 * cmH2O2mmHg /1000 # Mouth-Larynx resistance (mmHg.s/ml: PRU)
Cl = 1.27 / cmH2O2mmHg # Larynx compliance (ml/mmHg)
V₀l = 34.4 # Larynx zero pressure volume (ml)
Rlt = 0.3369 * cmH2O2mmHg /1000 # Larynx-Trachea resistance (mmHg.s/ml: PRU)
Ctr = 2.38 / cmH2O2mmHg # Trachea compliance (ml/mmHg)
V₀tr = 6.63 # Trachea zero pressure volume (ml)
Rtb = 0.3063 * cmH2O2mmHg /1000 # Trachea-Bronchial resistance (mmHg.s/ml: PRU)
Cb = 13.1 / cmH2O2mmHg # Bronchial compliance (ml/mmHg)
V₀b = 18.7 # Brochea zero pressure volume (ml)
RbA = 0.0817 * cmH2O2mmHg /1000 # Bronchea-Alveolar resistance (mmHg.s/ml: PRU)
CA = 200 / cmH2O2mmHg # Alveolar compliance (ml/mmHg)
Ccw = 244.5 / cmH2O2mmHg # Chest wall compliance (ml/mmHg)
FRC = 2400 # Functional Residual Capacity (ml)

V₀A = FRC + (CA * pₚₗ₀) - (V₀l) - (Ctr * (p₀ - pₚₗ₀) + V₀tr) - (Cb * (p₀ - pₚₗ₀) + V₀b) # Alveolar zero pressure volume (ml)

"""
Tissue Gas Exchange Model
"""

#### Tissue Volumes
Vₜ_brain = 1300 # Brain tissue volume (ml)
Vₜ_heart = 284 # Heart tissue volume (ml)
Vₜ_splanchnic = 2673 # Lung tissue volume (ml)
Vₜ_renal = 262 # Renal tissue volume (ml)
Vₜ_sm = 31200 # Skeletal muscle tissue volume (ml)
ub_mass_frac = 0.429 # Fraction of skeletal muscle in the upper body (Jansen 2000, https://journals.physiology.org/doi/full/10.1152/jappl.2000.89.1.81)
Vₜ_ub = Vₜ_sm * ub_mass_frac # Upper body (skeletal muscle) tissue volume (ml)
Vₜ_legs = Vₜ_sm * (1-ub_mass_frac) # Leg (skeletal muscle) tissue volume (ml)

#### Tissue O₂ Consumption Rate
MO₂_total = 250 # Total O₂ consumption rate (ml/min)
MO₂_spes_rat = 7.384 # Splanchnic to extrasplanchnic O₂ consumption rate ratio (Albanese, 2016)

MO₂_brain = 47.502/60 # Brain O₂ consumption rate (ml/s)
MO₂_heart = 24/60 # Heart O₂ consumption rate (ml/s) (Dynamic)
MO₂_sm = 51.6/60 # Skeletal muscle O₂ consumption rate (ml/s)
MO₂_ub = MO₂_sm * ub_mass_frac # Upper body (skeletal muscle) O₂ consumption rate (ml/s)
MO₂_legs = MO₂_sm * (1-ub_mass_frac) # Leg (skeletal muscle) O₂ consumption rate (ml/s)
MO₂_rem = MO₂_total/60 - (MO₂_brain + MO₂_heart + MO₂_ub + MO₂_legs) # Remaining O₂ consumption rate (ml/min)
MO₂_renal = MO₂_rem / (1 + MO₂_spes_rat) # Renal O₂ consumption rate (ml/s)
MO₂_splanchnic = MO₂_renal * MO₂_spes_rat # Splanchnic O₂ consumption rate (ml/s)

#### Tissue CO₂ Production Rate
RQ₀ = 0.84 # Respiratory Quotient (CO₂ produced / O₂ consumed)
MCO₂_total = MO₂_total * RQ₀ # Total CO₂ production rate (ml/min)

#### Heart Workload
Whₙₒₘ = 12660 # Nominal heart workload (mmHg.ml/s)
τ_heart_w = 5 # Time constant for heart workload (s)

"""
Lung Gas Exchange Model
"""
#### Environmental Conditions
FIO₂ = 21.0379/100 # Fraction of O₂ in the inspired air
FICO₂ = 0.0421/100 # Fraction of CO₂ in the inspired air (%)
K = 1.2103 # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
pₐₜₘ = 760.0 # Atmospheric pressure (mmHg)
p_ws = 47.0 # Water vapor pressure (mmHg)

#### O₂ Dissociation Curve
CₛₐₜO₂ = 9.0 * 0.0227 # Saturation of O₂ in the lungs (mmol/l) -> mlO2/ml
h₁ = 0.3836
α₁ = 0.03198 # (/mmHg)
β₁ = 0.008275 # (/mmHg)
K₁ = 14.99 # (mmHg)

#### CO₂ Dissociation Curve
CₛₐₜCO₂ = 86.11 *0.0227 # Saturation of CO₂ in the lungs (mmol/l) -> mlCO2/ml
h₂ = 1.819
α₂ = 0.05591 # (/mmHg)
β₂ = 0.03255 # (/mmHg)
K₂ = 194.4 # (mmHg)

#### Physiological Status
Hgb = 15.0 # Hemoglobin concentration (g/dl)
sh = 0.017 # Pulmonary shunt fraction (fraction of blood that bypasses the lungs)
DO₂ = 0.00139 # O₂ diffusion coefficient (cm^2/s)
DCO₂ = 0.00139 # CO₂ diffusion coefficient (cm^2/s)

#### Constants
Hgb_O₂_binding = 1.34 # O₂ capacity (ml O₂/g Hgb), 1.34
sol_O₂ = 0.003/100 # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)

DmO₂ = 22.7 * 22.7/60*1000/Pa2mmHg# Diffusion coefficient of O₂ in blood (mmol/min/kPa) -> ml/s/mmHg
DmCO₂ =  20 * DmO₂# CO₂ diffusion coefficient ml/s/mmHg
DlCO₂ = DmCO₂
DlO₂ = DmO₂



"""
Local Autoregulation (Magosso 2001)
"""

PaCO₂n = 40.0 # Set point for paCO₂ (mmHg)

Gbpn = 0.15 # ml/(mmHg.s)

CvbO₂n = 0.14
gbO₂ = 140 # Raised from 10 by Albanese, ml blood/mlO₂
A_auto = 20.9
B_auto = 92.8
C_auto = 10570
D_auto = -5.251

#Rhpn = 19.71 # mmHg.s/ml, same as Rcc
CvhO₂n = 0.11
ghO₂ = 490 # Raised from 35 by Albanese, ml blood/mlO₂
khCO₂ = 11.11 # mmHg

# Rmpn = 4.48 # mmHg.s/ml, replaced by leg and arm microvascular resistances
CvmO₂n = 0.155
gmO₂ = 420 # Raised from 30 by Albanese, ml blood/mlO₂
kmCO₂ = 142.86 # mmHg

τO₂ = 10 # s
τCO₂ = 20 # s

"""
BaroReflex Parameters
"""

τzb = 6.37 # (s)
τpb = 2.076 # (s)
fabₘᵢₙ = 2.52 # (spikes/s)
fabₘₐₓ = 47.78 # (spikes/s)
Pn = 92 # (mmHg)
kab = 11.76 # (mmHg)

"""
Cardiopulmonary Reflex Parameters
"""

τzr = 6.37
τpr = 2.076
Prn = 6
kcpr = 11.76/18*5
fcprₘₐₓ = (47.78 - 2.52) * 5 / 18
fcprₘᵢₙ = 2.52

"""
Peripheral Chemoreceptors (Ursino 2002)
"""

Ap = 600
Bp = 10.18
KO₂ = 200
KCO₂ = 1 # (/s)
Cₜ = 0.36 # (ml/ml)
Kstat = 20 # (/s)
τ_ph = 3.5 # (s)
τ_zh = 600 # (s)
τ_pl = 3.5 # (s)
Kdyn = 45 # (/s)

"""
Lung Stretch Receptors (Albanese 2016)
"""

τasr = 2 # (s)
Gasr = 0.01176 # (spikes/ml/s), modified from 23.29

"""
Pulmonary: Central Chemoreceptors (Albanese 2016)
"""

Dc = 8.0 # Central chemoreceptor delay (s)
G_cA = 0.850 * cmH2O2mmHg # Central chemoreceptor amplitude gain (mmHg/mmHg)
G_cf = 0.9 # Central chemoreceptor frequency gain (breaths/min/mmHg)
τ_cA = 105 # Amplitude time constant for central chemoreceptors (s)
τ_cf = 400 # Frequency time constant for central chemoreceptors (s)

"""
Pulmonary: Peripheral Chemoreceptors (Albanese 2016)
"""

Dp = 7.0 # Peripheral chemoreceptor delay (s)
G_pA = 1.310 * cmH2O2mmHg # Peripheral chemoreceptor amplitude gain (mmHg/ν)
G_pf = 0.8735 # Peripheral chemoreceptor frequency gain (breaths/min/ν)
fapc_set = 3.7 # Set point for fapc (spikes/s)
τ_pA = 83 # Amplitude time constant for peripheral chemoreceptors (s)
τ_pf = 147.78 # Frequency time constant for peripheral chemoreceptors (s)

"""
Unsorted Parameters
"""


# #### Afferent Chemoreflex
# τ_c = 2.0 # (s)
# facₘᵢₙ = 1.16 # (spikes/s)
# facₘₐₓ = 17.07 # (spikes/s)
# PO₂n = 45 # (mmHg)
# kac = 29.27 # (mmHg)


#### Efferent sympathetic pathway
fes∞ =2.1 # (spikes/s)
fes₀ = 16.11 # (spikes/s)
# fesₘᵢₙ = 2.66 # (spikes/s)
fesₘₐₓ = 60
kes = 0.0675 # (s)

Wbₛₚ = -1.1375 # Baroreflex to arterioles
Wcₛₚ = 1.716 # Chemoreceptors to arterioles
Wpₛₚ = -0.3997 # Pulmonary stretch receptors to arterioles

Wbₛᵥ = -1.0806 # Baroreflex to venous tone
Wcₛᵥ = 1.716 # Chemoreceptors to venous tone
Wpₛᵥ = -0.2907 # Pulmonary stretch receptors to venous tone


Wbₛₕ = -1.75 # Baroreflex to heart
Wcₛₕ = 1 # Modified from 5, Chemoreceptors to heart
Wpₛₕ = 0 # Pulmonary stretch receptors to heart




# New ???



#### Efferent vagal pathway
fev∞ = 6.3 # (spikes/s)
fev₀ = 3.2 # (spikes/s)
fab₀ = 25 # (spikes/s)
kev = 7.06 # (s)

Wcᵥ = 0.2
Wpᵥ = 0.103
θᵥ = -0.68 # (spikes/s)

#### CNS Hypoxic Response
χₘₐₓₛₚ = 13.32 # (spikes/s)
χₘₐₓₛₕ = 3.59 # (spikes/s)
τisc = 30 # (s)
χₘᵢₙₛₚ = 7.33 # (spikes/s)
χₘᵢₙₛₕ = -49.38 # (spikes/s)
PO₂nₛₚ = 30 # (mmHg)
PO₂nₛₕ = 45 # (mmHg)
kiscₛₚ = 2 # (mmHg)
kiscₛₕ = 6 # (mmHg)

#### Ventilatory Response
Gv = 0.125 # (l/ν)
τv = 3 # (s)
Dv = 6 # (s)
VTₙ = 0.583 # (l)
facₙ = 3.6 # (spikes/s)

#### Reflex Effectors
GEmaxlv = 0.475 # (mmHg/ml/ν)
GEmaxrv = 0.282 # (mmHg/ml/ν)
GRsp = 0.695 # (mmHg.s/ml/ν)
GRep = 0.695 # (mmHg.s/ml/ν)
GRmp = 0.695 # (mmHg.s/ml/ν)
Gvusv = -265.4 # (ml/ν)
Gvuev = -74.21 # (ml/ν)
Gvumv = -58.29 # (ml/ν)
GTs = -0.13 # (s/ν)
GTv = 0.09 # (s/ν)

τEmaxlv = 8 # (s)
τEmaxrv = 8 # (s)
τRsp = 6 # (s)
τRep = 6 # (s)
τRmp = 6 # (s)
τVusv = 20 # (s)
τVuev = 20 # (s)
τVumv = 20 # (s)
τTs = 2 # (s)
τTv = 1.5 # (s)

DEmaxlv = 2 # (s)
DEmaxrv = 2 # (s)
DRsp = 2 # (s)
DRep = 2 # (s)
DRmp = 2 # (s)
DVusv = 5 # (s)
DVuev = 5 # (s)
DVumv = 5 # (s)
DTs = 2 # (s)
DTv = 0.2 # (s)

Emaxlv₀ = 2.392 # (mmHg/ml)
Emaxrv₀ = 1.412 # (mmHg/ml)
Rsp₀ = 2.49 # (mmHg.s/ml)
Rep₀ = 1.655 # (mmHg.s/ml)
Rmp₀ = 2.106 # (mmHg.s/ml)
Vusv₀ = 1435.4 # (ml)
Vuev₀ = 640.73 # (ml)
Vumv₀ = 503.26 # (ml)
T₀ = 0.58 # (s)







for name in names(@__MODULE__; all=true, imported=false)
  # Only export if it's not a function or macro and is defined
  if isdefined(@__MODULE__, name) && !(name in (:eval, :include, :__doc__)) && !(name in names(Base))
      obj = getfield(@__MODULE__, name)
      # Only export if it's a variable, not a function
      if !(obj isa Function) && !(obj isa Type)
          @eval export $name
      end
  end
end

end