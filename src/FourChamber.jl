module FourChamber

#########################
# Preamble
#########################

using OrdinaryDiffEq
using DifferentialEquations
using ModelingToolkit
using Plots
using CSV
using DataFrames

#########################
# Element and Compartment Definitions
#########################
include("CompartmentDefinitions.jl")

#########################
# Parameters
#########################
include("ModelParameters.jl")

########################
# Instance Compartments
#########################
# Heart
@named LA = HeldtChamber(V₀=v0_la, p₀=p0_la, Eₘᵢₙ=Ed_la, Eₘₐₓ=Ees_la, τ=τ, τₑₛ=τes_la, Eshift=τav)
@named LV = HeldtChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Ed_lv, Eₘₐₓ=Ees_lv, τ=τ, τₑₛ=τes_lv, Eshift=0.0)
@named RA = HeldtChamber(V₀=v0_ra, p₀=p0_ra, Eₘᵢₙ=Ed_ra, Eₘₐₓ=Ees_ra, τ=τ, τₑₛ=τes_ra, Eshift=τav)
@named RV = HeldtChamber(V₀=v0_rv, p₀=p0_rv, Eₘᵢₙ=Ed_rv, Eₘₐₓ=Ees_rv, τ=τ, τₑₛ=τes_rv, Eshift=0.0)
@named MV = ResistorDiode(R=R_mv)
@named AV = ResistorDiode(R=R_av)
@named TV = ResistorDiode(R=R_tv)
@named PV = ResistorDiode(R=R_pv)
# Pulmonary Compartments
@named PulmArt = CRL(L.L=Lpa, R.R=Rpa, C.C=Cpa, C.p₀=pₜₕ, C.V₀=v0pa)
@named PulmCap = Resistor(R=Rpc)
@named PulmVein = StandardComp(R.R=Rpv, C.C=Cpv, C.p₀=pₜₕ, C.V₀=v0pv)
#Systemic Compartments
@named Body = StandardComp(R.R=Rsys, C.C=Csys, C.p₀=p₀, C.V₀=v0sys)

#########################
# Define Model
#########################
circ_eqs = [
  connect(LA.out, MV.in)
  connect(MV.out, LV.in)
  connect(LV.out, AV.in)
  connect(AV.out, Body.in)
  connect(Body.out, RA.in)
  connect(RA.out, TV.in)
  connect(TV.out, RV.in)
  connect(RV.out, PV.in)
  connect(PV.out, PulmArt.in)
  connect(PulmArt.out, PulmCap.in)
  connect(PulmCap.out, PulmVein.in)
  connect(PulmVein.out, LA.in)
]

#########################
# Define ODEs
#########################
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [LA, LV, RA, RV, MV, AV, TV, PV, PulmArt, PulmCap, PulmVein, Body])
circ_sys = structural_simplify(circ_model)
equations(expand_connections(circ_sys))
unknowns(circ_sys)

u0 = [
  LA.V => la0v
  LV.V => lv0v
  RA.V => ra0v
  RV.V => rv0v
  PulmArt.C.p => pt0pa
  PulmArt.L.q => qt0pa
  PulmVein.C.p => pt0pv
  Body.C.p => pt0sys
]


#########################
# Solve and Plot
#########################
prob = ODEProblem(circ_sys, u0, tspan)
Sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=17:0.01:20)
plot(Sol, idxs=[circ_model.LV.p, circ_model.LA.p, circ_model.PulmArt.C.p])

end # module OneChamber
