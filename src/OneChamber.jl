module OneChamber

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
@named LV = HeldtChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Ed_lv, Eₘₐₓ=Ees_lv, τ=τ, τₑₛ=τes_lv, Eshift=0.0)
@named R_aortic = ResistorDiode(R=R_av)
@named Body = StandardComp(R.R=Rsys, C.C=Csys)

#########################
# Define Model
#########################
circ_eqs = [
  connect(LV.out, R_aortic.in)
  connect(R_aortic.out, Body.in)
  connect(Body.out, LV.in)
]

#########################
# Define ODEs
#########################
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model, [LV, R_aortic, Body])
circ_sys = structural_simplify(circ_model)
equations(expand(circ_sys))
unknowns(circ_sys)

u0 = [
  Body.C.p => pt0sys
  LV.V => lv0v
]


#########################
# Solve and Plot
#########################
prob = ODEProblem(circ_sys, u0, tspan)
Sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=17:0.01:20)
plot(Sol, idxs=[circ_model.LV.p, circ_model.Body.C.p, circ_model.LV.V])

end # module OneChamber
