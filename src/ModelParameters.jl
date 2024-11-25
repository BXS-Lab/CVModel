#########################
# Global Parameters
#########################
@independent_variables t
tspan = (0, 20)


HRₙₒₘ = 60.0
pₜₕ = -4.0 # Thoracic Pressure
p₀ = 0.0  # External Pressure

τ = 1.0
Eshift = 0.0
Ev = Inf

#########################
# Cardiac Parameters
#########################
RR₀ = 60.0 / HRₙₒₘ
HR = HRₙₒₘ

#########################
# Cardiac Parameters
#########################

RR = 60.0 / HR
τav = 0.12 * sqrt(RR)
#### LA chamber parameters
v0_la = 24.0
p0_la = pₜₕ
Ed_la = 0.50
Ees_la = 0.61
τes_la = 0.2 * sqrt(RR)
#### LV chamber parameters
v0_lv = 55.0
p0_lv = pₜₕ
Ed_lv = 0.13
Ees_lv = 2.5
τes_lv = 0.3 * sqrt(RR)
#### RA chamber parameters
v0_ra = 14.0
p0_ra = pₜₕ
Ed_ra = 0.30
Ees_ra = 0.74
τes_ra = 0.2 * sqrt(RR)
#### RV chamber parameters
v0_rv = 46.0
p0_rv = pₜₕ
Ed_rv = 0.07
Ees_rv = 1.3
τes_rv = 0.3 * sqrt(RR)
## Mitral Valve parameters
R_mv = 0.01
## Aortic Valve parameters
R_av = 0.1
## Tricuspid Valve parameters
R_tv = 0.006
## Pulmonary Valve parameters
R_pv = 0.1

#########################
# Pulmonary Parameters
#########################
## Pulmonary Artery parameters
Rpa = 0.006
Cpa = 3.4
v0pa = 160.0
Lpa = 0.0017
## Pulmonary Capillary parameters
Rpc = 0.07
## Pulmonary Vein parameters
Rpv = 0.006
Cpv = 9.0
v0pv = 430.0

#########################
# Systemic Parameters
#########################

## Systemic parameters
Csys = 22.0
Rsys = 1.0
v0sys = 4200.0

#########################
# Initial Conditions
#########################

## Initial conditions
la0v = 24
lv0v = 55.0
ra0v = 14.0
rv0v = 46.0
pt0pa = 15.0
qt0pa = 0.0
pt0pv = 5.0
pt0sys = 90.0