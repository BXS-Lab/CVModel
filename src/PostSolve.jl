"""
This file contains the post-processing code for the cardiovascular model. Here some lumped values are calculated along with the beat-to-beat values for important systemic parameters. The beat-averaged compartment volumes are also calculated.

BXS Lab, UC Davis. Last updated May 2nd, 2025.
"""

"""
Calculate Lumped Elements
This code section calculates some of the lumped values including the total vascular volume, the instantaneous zero-pressure filling volume, the interstitial volume, the total volume (use as a check on conservation), and the total peripheral resistance (TPR). It also calculates the instantaneous heart rate (HR) based on the RR interval.
"""

Vvessel = Asc_A.C.V + BC_A.C.V + UpBd_art.C.V + UpBd_vein.C.V + SVC.C.V + Thor_A.C.V + Abd_A.C.V + Renal_art.C.V + Renal_vein.C.V + Splanchnic_art.C.V + Splanchnic_vein.C.V + Leg_art.C.V + Leg_vein.C.V + Abd_veins.C.V + Thor_IVC.C.V + RA.V + RV.V + Pulm_art.C.V + Pulm_vein.C.V + LA.V + LV.V

Vzpf = v0_Asc_A + v0_BC_A + v0_UpBd_art + UpBd_vein.C.V₀eff + v0_SVC + v0_Thor_A + v0_Abd_A + v0_Renal_art + Renal_vein.C.V₀eff + v0_Splanchnic_art + Splanchnic_vein.C.V₀eff + v0_Leg_art + Leg_vein.C.V₀eff + v0_Abd_veins + v0_Thor_IVC + v0_ra + v0_rv + v0pa + v0pv + v0_la + v0_lv

Vinterstitial = Interstitial.Vint
Vtotal = Vvessel + Vinterstitial
TPR = 1 / (1/UpBd_cap.R + 1/Renal_cap.R + 1/Splanchnic_cap.R + 1/Leg_cap.R)
HR = 60 / SA.RR_held
RR = SA.RR_held * 1000 # Convert to ms

"""
Calculate Beat-to-Beat values
This code section calculates the beat-tob-beat values for various model parameters:
SBP, DBP, MAP_integrated, MAP_formula, PP, CVP, CO, SV, SV_integrated, EF.
Plot using beat_times on the x-axis.
"""

ϕ_vals = Sol[SA.ϕ_wrapped] # Cardiac cycle phase signal
t_vals = Sol.t
Arterial_Pressure = Sol[Asc_A.pₜₘ]
LV_Volume = Sol[LV.V]
RA_Ptm = Sol[RA.pₜₘ]
Qlvo = -Sol[LV.out.q]
HR_beat = Sol[HR]
TPR_beat = Sol[TPR]

beat_indices = findall(diff(ϕ_vals) .< -0.5)  # reset points where ϕ jumps from ~1 back to ~0
push!(beat_indices, length(t_vals))  # include end of last beat

beat_times = Float64[]
SBP = Float64[]
DBP = Float64[]
MAP_integrated = Float64[]
CVP = Float64[]
Qlvo_mean = Float64[]
CO = Float64[]
LV_Volume_max = Float64[]
LV_Volume_min = Float64[]
SV_integrated = Float64[]
HR_beats = Float64[]
TPR_beats = Float64[]

function compute_beat_metrics(beat_indices, t_vals, Arterial_Pressure, LV_Volume, RA_Ptm, Qlvo, HR_beat)
  start_idx = 1
  for stop_idx in beat_indices
      beat_range = start_idx:stop_idx
      push!(SBP, maximum(Arterial_Pressure[beat_range]))
      push!(DBP, minimum(Arterial_Pressure[beat_range]))
      push!(LV_Volume_max, maximum(LV_Volume[beat_range]))
      push!(LV_Volume_min, minimum(LV_Volume[beat_range]))
      push!(beat_times, t_vals[beat_range][1])
      push!(HR_beats, HR_beat[beat_range][1])
      push!(TPR_beats, TPR_beat[beat_range][1])
      t_segment = t_vals[beat_range]
      p_segment = Arterial_Pressure[beat_range]
      ra_segment = RA_Ptm[beat_range]
      co_segment = Qlvo[beat_range]
      dt = diff(t_segment)
      integrand_map = p_segment[1:end-1] .* dt
      integrand_ra = ra_segment[1:end-1] .* dt
      integrand_co = co_segment[1:end-1] .* dt
      push!(MAP_integrated, sum(integrand_map) / (t_segment[end] - t_segment[1]))
      push!(CVP, sum(integrand_ra) / (t_segment[end] - t_segment[1]))
      push!(Qlvo_mean, sum(integrand_co) / (t_segment[end] - t_segment[1]))
      push!(SV_integrated, sum(integrand_co))
      start_idx = stop_idx + 1
  end
end
compute_beat_metrics(beat_indices, t_vals, Arterial_Pressure, LV_Volume, RA_Ptm, Qlvo, HR_beat)

PP = SBP .- DBP
MAP_formula = DBP .+ (1/3) .* PP
CO = Qlvo_mean .* 6/100  # Convert to L/min
SV = LV_Volume_max .- LV_Volume_min
EF = SV ./ LV_Volume_max  # Ejection fraction = Stroke volume / EDV
CO_Calculated = SV .* HR_beats ./ 1000
RR_beats = 60 ./ HR_beats *1000

"""
Calculate Beat-Averaged Compartment Volumes
This code section calculates the average compartment volumes over the cardiac cycle.
"""

Asc_A_V = Sol[Asc_A.C.V]
BC_A_V = Sol[BC_A.C.V]
UpBd_art_V = Sol[UpBd_art.C.V]
UpBd_vein_V = Sol[UpBd_vein.C.V]
SVC_V = Sol[SVC.C.V]
Thor_A_V = Sol[Thor_A.C.V]
Abd_A_V = Sol[Abd_A.C.V]
Renal_art_V = Sol[Renal_art.C.V]
Renal_vein_V = Sol[Renal_vein.C.V]
Splanchnic_art_V = Sol[Splanchnic_art.C.V]
Splanchnic_vein_V = Sol[Splanchnic_vein.C.V]
Leg_art_V = Sol[Leg_art.C.V]
Leg_vein_V = Sol[Leg_vein.C.V]
Abd_veins_V = Sol[Abd_veins.C.V]
Thor_IVC_V = Sol[Thor_IVC.C.V]
RA_V = Sol[RA.V]
RV_V = Sol[RV.V]
Pulm_art_V = Sol[Pulm_art.C.V]
Pulm_vein_V = Sol[Pulm_vein.C.V]
LA_V = Sol[LA.V]
LV_V = Sol[LV.V]

Asc_A_Vmean = Float64[]
BC_A_Vmean = Float64[]
UpBd_art_Vmean = Float64[]
UpBd_vein_Vmean = Float64[]
SVC_Vmean = Float64[]
Thor_A_Vmean = Float64[]
Abd_A_Vmean = Float64[]
Renal_art_Vmean = Float64[]
Renal_vein_Vmean = Float64[]
Splanchnic_art_Vmean = Float64[]
Splanchnic_vein_Vmean = Float64[]
Leg_art_Vmean = Float64[]
Leg_vein_Vmean = Float64[]
Abd_veins_Vmean = Float64[]
Thor_IVC_Vmean = Float64[]
RA_Vmean = Float64[]
RV_Vmean = Float64[]
Pulm_art_Vmean = Float64[]
Pulm_vein_Vmean = Float64[]
LA_Vmean = Float64[]
LV_Vmean = Float64[]

function compute_beat_averaged_volumes(beat_indices, t_vals, Asc_A_V, BC_A_V, UpBd_art_V, UpBd_vein_V, SVC_V, Thor_A_V, Abd_A_V, Renal_art_V, Renal_vein_V, Splanchnic_art_V, Splanchnic_vein_V, Leg_art_V, Leg_vein_V, Abd_veins_V, Thor_IVC_V, RA_V, RV_V, Pulm_art_V, Pulm_vein_V, LA_V, LV_V)
  start_idx = 1
  for stop_idx in beat_indices
      beat_range = start_idx:stop_idx
      t_segment = t_vals[beat_range]
      Asc_A_v_segment = Asc_A_V[beat_range]
      BC_A_v_segment = BC_A_V[beat_range]
      UpBd_art_v_segment = UpBd_art_V[beat_range]
      UpBd_vein_v_segment = UpBd_vein_V[beat_range]
      SVC_v_segment = SVC_V[beat_range]
      Thor_A_v_segment = Thor_A_V[beat_range]
      Abd_A_v_segment = Abd_A_V[beat_range]
      Renal_art_v_segment = Renal_art_V[beat_range]
      Renal_vein_v_segment = Renal_vein_V[beat_range]
      Splanchnic_art_v_segment = Splanchnic_art_V[beat_range]
      Splanchnic_vein_v_segment = Splanchnic_vein_V[beat_range]
      Leg_art_v_segment = Leg_art_V[beat_range]
      Leg_vein_v_segment = Leg_vein_V[beat_range]
      Abd_veins_v_segment = Abd_veins_V[beat_range]
      Thor_IVC_v_segment = Thor_IVC_V[beat_range]
      RA_v_segment = RA_V[beat_range]
      RV_v_segment = RV_V[beat_range]
      Pulm_art_v_segment = Pulm_art_V[beat_range]
      Pulm_vein_v_segment = Pulm_vein_V[beat_range]
      LA_v_segment = LA_V[beat_range]
      LV_v_segment = LV_V[beat_range]
      dt = diff(t_segment)
      integrand_Asc_A = Asc_A_v_segment[1:end-1] .* dt
      integrand_BC_A = BC_A_v_segment[1:end-1] .* dt
      integrand_UpBd_art = UpBd_art_v_segment[1:end-1] .* dt
      integrand_UpBd_vein = UpBd_vein_v_segment[1:end-1] .* dt
      integrand_SVC = SVC_v_segment[1:end-1] .* dt
      integrand_Thor_A = Thor_A_v_segment[1:end-1] .* dt
      integrand_Abd_A = Abd_A_v_segment[1:end-1] .* dt
      integrand_Renal_art = Renal_art_v_segment[1:end-1] .* dt
      integrand_Renal_vein = Renal_vein_v_segment[1:end-1] .* dt
      integrand_Splanchnic_art = Splanchnic_art_v_segment[1:end-1] .* dt
      integrand_Splanchnic_vein = Splanchnic_vein_v_segment[1:end-1] .* dt
      integrand_Leg_art = Leg_art_v_segment[1:end-1] .* dt
      integrand_Leg_vein = Leg_vein_v_segment[1:end-1] .* dt
      integrand_Abd_veins = Abd_veins_v_segment[1:end-1] .* dt
      integrand_Thor_IVC = Thor_IVC_v_segment[1:end-1] .* dt
      integrand_RA = RA_v_segment[1:end-1] .* dt
      integrand_RV = RV_v_segment[1:end-1] .* dt
      integrand_Pulm_art = Pulm_art_v_segment[1:end-1] .* dt
      integrand_Pulm_vein = Pulm_vein_v_segment[1:end-1] .* dt
      integrand_LA = LA_v_segment[1:end-1] .* dt
      integrand_LV = LV_v_segment[1:end-1] .* dt
      push!(Asc_A_Vmean, sum(integrand_Asc_A) / (t_segment[end] - t_segment[1]))
      push!(BC_A_Vmean, sum(integrand_BC_A) / (t_segment[end] - t_segment[1]))
      push!(UpBd_art_Vmean, sum(integrand_UpBd_art) / (t_segment[end] - t_segment[1]))
      push!(UpBd_vein_Vmean, sum(integrand_UpBd_vein) / (t_segment[end] - t_segment[1]))
      push!(SVC_Vmean, sum(integrand_SVC) / (t_segment[end] - t_segment[1]))
      push!(Thor_A_Vmean, sum(integrand_Thor_A) / (t_segment[end] - t_segment[1]))
      push!(Abd_A_Vmean, sum(integrand_Abd_A) / (t_segment[end] - t_segment[1]))
      push!(Renal_art_Vmean, sum(integrand_Renal_art) / (t_segment[end] - t_segment[1]))
      push!(Renal_vein_Vmean, sum(integrand_Renal_vein) / (t_segment[end] - t_segment[1]))
      push!(Splanchnic_art_Vmean, sum(integrand_Splanchnic_art) / (t_segment[end] - t_segment[1]))
      push!(Splanchnic_vein_Vmean, sum(integrand_Splanchnic_vein) / (t_segment[end] - t_segment[1]))
      push!(Leg_art_Vmean, sum(integrand_Leg_art) / (t_segment[end] - t_segment[1]))
      push!(Leg_vein_Vmean, sum(integrand_Leg_vein) / (t_segment[end] - t_segment[1]))
      push!(Abd_veins_Vmean, sum(integrand_Abd_veins) / (t_segment[end] - t_segment[1]))
      push!(Thor_IVC_Vmean, sum(integrand_Thor_IVC) / (t_segment[end] - t_segment[1]))
      push!(RA_Vmean, sum(integrand_RA) / (t_segment[end] - t_segment[1]))
      push!(RV_Vmean, sum(integrand_RV) / (t_segment[end] - t_segment[1]))
      push!(Pulm_art_Vmean, sum(integrand_Pulm_art) / (t_segment[end] - t_segment[1]))
      push!(Pulm_vein_Vmean, sum(integrand_Pulm_vein) / (t_segment[end] - t_segment[1]))
      push!(LA_Vmean, sum(integrand_LA) / (t_segment[end] - t_segment[1]))
      push!(LV_Vmean, sum(integrand_LV) / (t_segment[end] - t_segment[1]))
      start_idx = stop_idx + 1
  end
end
compute_beat_averaged_volumes(beat_indices, t_vals, Asc_A_V, BC_A_V, UpBd_art_V, UpBd_vein_V, SVC_V, Thor_A_V, Abd_A_V, Renal_art_V, Renal_vein_V, Splanchnic_art_V, Splanchnic_vein_V, Leg_art_V, Leg_vein_V, Abd_veins_V, Thor_IVC_V, RA_V, RV_V, Pulm_art_V, Pulm_vein_V, LA_V, LV_V)