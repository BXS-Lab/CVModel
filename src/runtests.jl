##
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq
using Test
using CSV
using DataFrames
##
# @testset "WK5" begin

# end


    ##
    include("ShiParam.jl")

    @independent_variables t

    ## Ventricles
    @named LV = ShiChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    @named LA = ShiChamber(V₀=v0_la, p₀=p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la / 2, τₑₚ=τpww_la, Eshift=τpwb_la)
    @named RV = ShiChamber(V₀=v0_rv, p₀=p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    # @named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)
    @named RA = ShiAtrium(V₀=v0_ra, p₀=1, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τpwb=τpwb_ra, τpww=τpww_ra) #, Ev=Inf)

    ## Valves as simple valves
    @named AV = OrificeValve(CQ=CQ_AV)
    @named MV = OrificeValve(CQ=CQ_MV)
    @named TV = OrificeValve(CQ=CQ_TV)
    @named PV = OrificeValve(CQ=CQ_PV)

    ####### Systemic Loop #######
    # Systemic Aortic Sinus ##
    @named SAS = CRL(C=Csas, R=Rsas, L=Lsas)
    # Systemic Artery ##
    @named SAT = CRL(C=Csat, R=Rsat, L=Lsat)
    # Systemic Arteriole ##
    @named SAR = Resistor(R=Rsar)
    # Systemic Capillary ##
    @named SCP = Resistor(R=Rscp)
    # Systemic Vein ##
    @named SVN = CR(R=Rsvn, C=Csvn)

    ####### Pulmonary Loop #######
    # Pulmonary Aortic Sinus ##
    @named PAS = CR(C=Cpas, R=Rpas)
    # Pulmonary Artery ##
    @named PAT = CR(C=Cpat, R=Rpat)
    # Pulmonary Arteriole ##
    @named PAR = Resistor(R=Rpar)
    # Pulmonary Capillary ##
    @named PCP = Resistor(R=Rpcp)
    # Pulmonary Vein ##
    @named PVN = CR(R=Rpvn, C=Cpvn)

    ##
    circ_eqs = [
        connect(LV.out, AV.in)
        connect(AV.out, SAS.in)
        connect(SAS.out, SAT.in)
        connect(SAT.out, SAR.in)
        connect(SAR.out, SCP.in)
        connect(SCP.out, SVN.in)
        connect(SVN.out, RA.in)
        connect(RA.out, TV.in)
        connect(TV.out, RV.in)
        connect(RV.out, PV.in)
        connect(PV.out, PAS.in)
        connect(PAS.out, PAT.in)
        connect(PAT.out, PAR.in)
        connect(PAR.out, PCP.in)
        connect(PCP.out, PVN.in)
        connect(PVN.out, LA.in)
        connect(LA.out, MV.in)
        connect(MV.out, LV.in)
    ]

    ## Compose the whole ODE system
    @named _circ_model = ODESystem(circ_eqs, t)
    @named circ_model = compose(_circ_model,
        [LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

    ## And simplify it
    circ_sys = structural_simplify(circ_model)

    unknowns(circ_sys)
    ## Setup ODE
    # Initial Conditions for Shi Valve
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0,0, 0, 0, 0]
    # and for OrificeValve --- Commment this next line to use ShiValves
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

    u0 = [
        LV.V => LV_Vt0
        RV.V => RV_Vt0
        LA.V => LA_Vt0
        RA.V => RA_Vt0
        SAS.C.p => pt0sas
        SAS.L.q => qt0sas
        SAT.C.p => pt0sat
        SAT.L.q => qt0sat
        SVN.C.p => pt0svn
        PAS.C.p => pt0pas
        PAS.L.q => qt0pas
        PAT.C.p => pt0pat
        PAT.L.q => qt0pat
        PVN.C.p => pt0pvn
    ]

    prob = ODEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time ShiSimpleSolV = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)
    # ShiSimpleSolV = ShiSimpleSolV(19:0.01:20)

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiSimple.csv", DataFrame)

    @test SciMLBase.successful_retcode(ShiSimpleSolV)
    @test sum((ShiSimpleSolV[LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
end

#