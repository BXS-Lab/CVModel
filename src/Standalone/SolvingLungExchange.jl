"""
Preamble
"""

using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using Symbolics
using NLsolve

@parameters t
D = Differential(t)

"""
Base Lung Model
"""

@mtkmodel LungGasExchange begin
  @parameters begin
    _FIO₂ = FIO₂ # Fractional concentration of O₂ in the inspired air (%)
    _FICO₂ = FICO₂ # Fractional concentration of CO₂ in the inspired air (%)
    _sh = sh # Shunt fraction
    _v0pa = v0pa # Unstressed volume of the pulmonary arterial blood (ml)
    _K = K # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
    _CₛₐₜO₂ = CₛₐₜO₂ # O₂ saturation constant (ml/ml)
    _CₛₐₜCO₂ = CₛₐₜCO₂ # CO₂ saturation constant (ml/ml)
    _h₁ = h₁ # O₂ saturation exponent
    _h₂ = h₂ # CO₂ saturation exponent
    _K₁ = K₁ # O₂ saturation constant (mmHg)
    _K₂ = K₂ # CO₂ saturation constant (mmHg)
    _α₁ = α₁ # O₂ saturation constant (/mmHg)
    _α₂ = α₂ # CO₂ saturation constant (/mmHg)
    _β₁ = β₁ # O₂ saturation constant (/mmHg)
    _β₂ = β₂ # CO₂ saturation constant (/mmHg)
    _pₐₜₘ = pₐₜₘ # Absolute Atmospheric pressure (mmHg)
    _p_ws = p_ws # Water vapor pressure (mmHg)
    _sol_O₂ = sol_O₂ # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)
    _Hgb_O₂_binding = Hgb_O₂_binding # Hemoglobin O₂ binding constant (ml O₂/ml blood/mmHg)
    _Hgb = Hgb # Hemoglobin concentration (g/dl)
    _DlO₂ = DlO₂ # Diffusion coefficient of O₂ in blood (ml O₂/ml blood/mmHg)
    _DlCO₂ = DlCO₂ # Diffusion coefficient of CO₂ in blood (ml CO₂/ml blood/mmHg)
  end
  @variables begin
    #### Connected from Lung Model
    Vrᵢₙ(t) # Flow rate (ml/s), must be connected externally
    Vr_A(t) # Alveolar flow rate (ml/s), must be connected externally
    V_D(t) # Dead space volume (ml), must be connected externally
    V_A(t) # Alveolar volume (ml), must be connected externally

    #### Connected from Pulmonary Artery
    qpa(t) # Pulmonary blood flow (ml/s), must be connected externally
    Vpa(t) # Pulmonary arterial blood volume (ml), must be connected externally
    cvO₂(t) # O₂ concentration in the deoxygenated blood (ml/ml), must be connected externally
    cvCO₂(t) # CO₂ concentration in the deoxygenated blood (ml/ml), must be connected externally


    FDO₂(t) # Fractional concentration of O₂ in the dead space (%)
    FDCO₂(t) # Fractional concentration of CO₂ in the dead space (%)
    FAO₂(t) # Fractional concentration of O₂ in the alveolar space (%)
    FACO₂(t) # Fractional concentration of CO₂ in the alveolar space (%)

    qpp(t) # Pulmonary blood flow (ml/s)
    qps(t) # Pulmonary shunted blood flow (ml/s)
    Vpp(t) # Pulmonary peripheral blood volume (ml)

    cppO₂(t) # O₂ concentration in the pulmonary blood (ml/ml)
    cppCO₂(t) # CO₂ concentration in the pulmonary blood (ml/ml)
    XppO₂(t) # O₂ partial pressure
    XppCO₂(t) # CO₂ partial pressure
    pppO₂(t) # O₂ partial pressure in pulmonary blood (mmHg)
    pppCO₂(t) # CO₂ partial pressure in pulmonary blood (mmHg)

    p_AO₂(t) # O₂ partial pressure in the alveolar space (mmHg)
    p_ACO₂(t) # CO₂ partial pressure in the alveolar space (mmHg)

    caO₂(t) # O₂ concentration in the arterial blood (ml/ml) (output)
    caCO₂(t) # CO₂ concentration in the arterial blood (ml/ml) (output)

    VrO₂p(t) # O₂ uptake (ml/s)
    VrCO₂p(t) # CO₂ uptake (ml/s)
    VrO₂d(t) # O₂ diffusion (ml/s)
    VrCO₂d(t) # CO₂ diffusion (ml/s)
    VrO₂(t) # O₂ uptake (ml/s)
    VrCO₂(t) # CO₂ uptake (ml/s)

    XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
    XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
    paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
    paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

    SaO₂(t) # O₂ saturation in the arterial blood (%)
  end
  @equations begin
    #### Conservation of mass equations
    D(FDO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FIO₂ - FDO₂), Vr_A * (FDO₂ - FAO₂)) / V_D
    D(FDCO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FICO₂ - FDCO₂), Vr_A * (FDCO₂ - FACO₂)) / V_D
    D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * D(cppO₂))) / V_A
    D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * D(cppCO₂))) / V_A
    Vpp ~ (Vpa - v0pa) * (1 - _sh) + v0pa


    # #### Dissociation equations
    cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
    XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
    cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
    XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

    # #### Instantaneous equilibrium equations
    p_AO₂ ~ pppO₂
    p_ACO₂ ~ pppCO₂

    #### Gas fraction to partial pressure relationships
    p_AO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
    p_ACO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

    #### Mixing between capillary and shunted blood
    caO₂ ~ (1 - _sh) * cppO₂ + _sh * cvO₂
    caCO₂ ~ (1 - _sh) * cppCO₂ + _sh * cvCO₂

    qpp ~ qpa * (1 - _sh)
    qps ~ qpa * _sh

    #### O₂ saturation in arterial blood
    XaO₂ ~ max(((caO₂ / _CₛₐₜO₂) / (1 - (caO₂ / _CₛₐₜO₂))),1e-8)^_h₁
    XaCO₂ ~ max(((caCO₂ / _CₛₐₜCO₂) / (1 - (caCO₂ / _CₛₐₜCO₂))),1e-8)^_h₂

    paO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) + _K₁*XaO₂*_β₂ - _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₂ + 2*_K₂*XaCO₂*_α₂*_β₁)

    paCO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) - _K₁*XaO₂*_β₂ + _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₁ + 2*_K₁*XaO₂*_α₁*_β₂)

    SaO₂ ~ (caO₂ - paO₂ * _sol_O₂) / (_Hgb * _Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
  end
end

"""
Define parameters
"""

ps = @parameters begin
    _FIO₂
    _FICO₂
    _sh
    _v0pa
    _K
    _CₛₐₜO₂
    _CₛₐₜCO₂
    _h₁
    _h₂
    _K₁
    _K₂
    _α₁
    _α₂
    _β₁
    _β₂
    _pₐₜₘ
    _p_ws
    _sol_O₂
    _Hgb_O₂_binding
    _Hgb
    _DlO₂
    _DlCO₂
end

sts = @variables begin
    Vrᵢₙ(t)
    Vr_A(t)
    V_D(t)
    V_A(t)

    qpa(t)
    Vpa(t)
    cvO₂(t)
    cvCO₂(t)

    FDO₂(t)
    FDCO₂(t)
    FAO₂(t)
    FACO₂(t)

    qpp(t)
    qps(t)
    Vpp(t)

    cppO₂(t)
    cppCO₂(t)
    XppO₂(t)
    XppCO₂(t)
    pppO₂(t)
    pppCO₂(t)

    p_AO₂(t)
    p_ACO₂(t)

    caO₂(t)
    caCO₂(t)

    VrO₂p(t)
    VrCO₂p(t)
    VrO₂d(t)
    VrCO₂d(t)
    VrO₂(t)
    VrCO₂(t)

    XaO₂(t)
    XaCO₂(t)
end

"""
Start with O₂
"""

#### Relevant Equations
D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * D(cppO₂))) / V_A
cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
pppO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
pppCO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

#### Find dpppO₂/dt and dpppCO₂/dt
D(pppO₂) ~ D(FAO₂) * (_pₐₜₘ - _p_ws)
D(pppCO₂) ~ D(FACO₂) * (_pₐₜₘ - _p_ws)

#### Find dXppO₂/dt
@variables pppO₂(t) pppCO₂(t)
@parameters _β₁ _K₁ _α₁

XppO₂_def = pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
d_XppO₂_dt = expand_derivatives(Differential(t)(XppO₂_def))

d_XppO₂_dt = (_β₁*pppO₂*D(pppCO₂) + (1 + _β₁*pppCO₂)*D(pppO₂)) / (_K₁*(1 + _α₁*pppCO₂)) - _K₁*_α₁*(((1 + _β₁*pppCO₂)*pppO₂) / ((_K₁^2)*((1 + _α₁*pppCO₂)^2)))*D(pppCO₂)

#### Write dXppO₂/dt in terms of FAO₂ and FACO₂
pppO₂ ~ (FAO₂ * (_pₐₜₘ - _p_ws))
pppCO₂ ~ (FACO₂ * (_pₐₜₘ - _p_ws))
D(pppO₂) ~ (D(FAO₂) * (_pₐₜₘ - _p_ws))
D(pppCO₂) ~ (D(FACO₂) * (_pₐₜₘ - _p_ws))

d_XppO₂_dt = (_β₁*(FAO₂ * (_pₐₜₘ - _p_ws))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + (1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₁*(1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))) - _K₁*_α₁*(((1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(FAO₂ * (_pₐₜₘ - _p_ws))) / ((_K₁^2)*((1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FACO₂) * (_pₐₜₘ - _p_ws))

#### Find dcppO₂/dXppO₂
@parameters XppO₂ _CₛₐₜO₂ _h₁

cppO₂_def = _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
d_cppO₂_d_XppO₂ = expand_derivatives(ModelingToolkit.derivative(cppO₂_def, XppO₂))

d_cppO₂_d_XppO₂ = (-(XppO₂^(-1 + 1 / _h₁))*((_CₛₐₜO₂*(XppO₂^(1 / _h₁))) / ((1 + XppO₂^(1 / _h₁))^2))) / _h₁ + (_CₛₐₜO₂*(XppO₂^(-1 + 1 / _h₁))) / (_h₁*(1 + XppO₂^(1 / _h₁)))

#### Write XppO₂ in terms of FAO₂ and FACO₂

XppO₂ ~ ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))

#### Find dcppO₂

d_cppO₂_d_t = ((-(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))*((_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))) / ((1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))^2))) / _h₁ + (_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))) / (_h₁*(1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁)))) * ((_β₁*(FAO₂ * (_pₐₜₘ - _p_ws))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + (1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₁*(1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))) - _K₁*_α₁*(((1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(FAO₂ * (_pₐₜₘ - _p_ws))) / ((_K₁^2)*((1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)))

#### Substitute into original equation for D(FAO₂)

D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * (((-(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))*((_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))) / ((1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))^2))) / _h₁ + (_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))) / (_h₁*(1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁)))) * ((_β₁*(FAO₂ * (_pₐₜₘ - _p_ws))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + (1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₁*(1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))) - _K₁*_α₁*(((1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(FAO₂ * (_pₐₜₘ - _p_ws))) / ((_K₁^2)*((1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)))))) / V_A

#### Rearrange to find D(FAO₂)

eq = D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * (((-(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))*((_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))) / ((1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁))^2))) / _h₁ + (_CₛₐₜO₂*(((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₁))) / (_h₁*(1 + ((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₁)))) * ((_β₁*(FAO₂ * (_pₐₜₘ - _p_ws))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + (1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₁*(1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))) - _K₁*_α₁*(((1 + _β₁*(FACO₂ * (_pₐₜₘ - _p_ws)))*(FAO₂ * (_pₐₜₘ - _p_ws))) / ((_K₁^2)*((1 + _α₁*(FACO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)))))) / V_A

solution = solve_for(eq, D(FAO₂))[1]

#### Final Result

D(FAO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-(cppO₂ - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁))

#### Really final results

D(FAO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁))

"""
Next do CO₂
"""

D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * D(cppCO₂))) / V_A
cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))
pppO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
pppCO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

#### Find dpppO₂/dt and dpppCO₂/dt
D(pppO₂) ~ D(FAO₂) * (_pₐₜₘ - _p_ws)
D(pppCO₂) ~ D(FACO₂) * (_pₐₜₘ - _p_ws)

#### Find dXppCO₂/dt
@variables pppO₂(t) pppCO₂(t)
@parameters _β₂ _K₂ _α₂

XppCO₂_def = pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))
d_XppCO₂_dt = expand_derivatives(Differential(t)(XppCO₂_def))

d_XppCO₂_dt = ((1 + _β₂*pppO₂)*D(pppCO₂) + _β₂*pppCO₂*D(pppO₂)) / (_K₂*(1 + _α₂*pppO₂)) - _K₂*_α₂*(((1 + _β₂*pppO₂)*pppCO₂) / ((_K₂^2)*((1 + _α₂*pppO₂)^2)))*D(pppO₂)


#### Write dXppCO₂/dt in terms of FAO₂ and FACO₂
pppO₂ ~ (FAO₂ * (_pₐₜₘ - _p_ws))
pppCO₂ ~ (FACO₂ * (_pₐₜₘ - _p_ws))
D(pppO₂) ~ (D(FAO₂) * (_pₐₜₘ - _p_ws))
D(pppCO₂) ~ (D(FACO₂) * (_pₐₜₘ - _p_ws))

d_XppCO₂_dt = ((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + _β₂*(FACO₂ * (_pₐₜₘ - _p_ws))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₂*(1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))) - _K₂*_α₂*(((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(FACO₂ * (_pₐₜₘ - _p_ws))) / ((_K₂^2)*((1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FAO₂) * (_pₐₜₘ - _p_ws))

#### Find dcppCO₂/dXppCO₂
@parameters XppCO₂ _CₛₐₜCO₂ _h₂

cppCO₂_def = _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
d_cppCO₂_d_XppCO₂ = expand_derivatives(ModelingToolkit.derivative(cppCO₂_def, XppCO₂))

d_cppCO₂_d_XppCO₂ = (-((_CₛₐₜCO₂*(XppCO₂^(1 / _h₂))) / ((1 + XppCO₂^(1 / _h₂))^2))*(XppCO₂^(-1 + 1 / _h₂))) / _h₂ + (_CₛₐₜCO₂*(XppCO₂^(-1 + 1 / _h₂))) / (_h₂*(1 + XppCO₂^(1 / _h₂)))

#### Write XppCO₂ in terms of FAO₂ and FACO₂

XppCO₂ ~ ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))

#### Find dcppCO₂

d_cppCO₂_d_t = ((-((_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))) / ((1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))^2))*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / _h₂ + (_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂)))) * (((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + _β₂*(FACO₂ * (_pₐₜₘ - _p_ws))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₂*(1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))) - _K₂*_α₂*(((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(FACO₂ * (_pₐₜₘ - _p_ws))) / ((_K₂^2)*((1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FAO₂) * (_pₐₜₘ - _p_ws)))

#### Substitute into original equation for D(FACO₂)

D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * (((-((_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))) / ((1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))^2))*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / _h₂ + (_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂)))) * (((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + _β₂*(FACO₂ * (_pₐₜₘ - _p_ws))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₂*(1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))) - _K₂*_α₂*(((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(FACO₂ * (_pₐₜₘ - _p_ws))) / ((_K₂^2)*((1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FAO₂) * (_pₐₜₘ - _p_ws)))))) / V_A

#### Rearrange to find D(FACO₂)

eq = D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * (((-((_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))) / ((1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂))^2))*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / _h₂ + (_CₛₐₜCO₂*(((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws)))))^(1 / _h₂)))) * (((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(D(FACO₂) * (_pₐₜₘ - _p_ws)) + _β₂*(FACO₂ * (_pₐₜₘ - _p_ws))*(D(FAO₂) * (_pₐₜₘ - _p_ws))) / (_K₂*(1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))) - _K₂*_α₂*(((1 + _β₂*(FAO₂ * (_pₐₜₘ - _p_ws)))*(FACO₂ * (_pₐₜₘ - _p_ws))) / ((_K₂^2)*((1 + _α₂*(FAO₂ * (_pₐₜₘ - _p_ws)))^2)))*(D(FAO₂) * (_pₐₜₘ - _p_ws)))))) / V_A

solution = solve_for(eq, D(FACO₂))[1]

#### Final Result

D(FACO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-(cppCO₂ - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂))

#### Really final results

D(FACO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂))



"""
Substitute Back into Original Lung Model
"""

@mtkmodel LungGasExchange begin
  @parameters begin
    _FIO₂ = FIO₂ # Fractional concentration of O₂ in the inspired air (%)
    _FICO₂ = FICO₂ # Fractional concentration of CO₂ in the inspired air (%)
    _sh = sh # Shunt fraction
    _v0pa = v0pa # Unstressed volume of the pulmonary arterial blood (ml)
    _K = K # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
    _CₛₐₜO₂ = CₛₐₜO₂ # O₂ saturation constant (ml/ml)
    _CₛₐₜCO₂ = CₛₐₜCO₂ # CO₂ saturation constant (ml/ml)
    _h₁ = h₁ # O₂ saturation exponent
    _h₂ = h₂ # CO₂ saturation exponent
    _K₁ = K₁ # O₂ saturation constant (mmHg)
    _K₂ = K₂ # CO₂ saturation constant (mmHg)
    _α₁ = α₁ # O₂ saturation constant (/mmHg)
    _α₂ = α₂ # CO₂ saturation constant (/mmHg)
    _β₁ = β₁ # O₂ saturation constant (/mmHg)
    _β₂ = β₂ # CO₂ saturation constant (/mmHg)
    _pₐₜₘ = pₐₜₘ # Absolute Atmospheric pressure (mmHg)
    _p_ws = p_ws # Water vapor pressure (mmHg)
    _sol_O₂ = sol_O₂ # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)
    _Hgb_O₂_binding = Hgb_O₂_binding # Hemoglobin O₂ binding constant (ml O₂/ml blood/mmHg)
    _Hgb = Hgb # Hemoglobin concentration (g/dl)
    _DlO₂ = DlO₂ # Diffusion coefficient of O₂ in blood (ml O₂/ml blood/mmHg)
    _DlCO₂ = DlCO₂ # Diffusion coefficient of CO₂ in blood (ml CO₂/ml blood/mmHg)
  end
  @variables begin
    #### Connected from Lung Model
    Vrᵢₙ(t) # Flow rate (ml/s), must be connected externally
    Vr_A(t) # Alveolar flow rate (ml/s), must be connected externally
    V_D(t) # Dead space volume (ml), must be connected externally
    V_A(t) # Alveolar volume (ml), must be connected externally

    #### Connected from Pulmonary Artery
    qpa(t) # Pulmonary blood flow (ml/s), must be connected externally
    Vpa(t) # Pulmonary arterial blood volume (ml), must be connected externally
    cvO₂(t) # O₂ concentration in the deoxygenated blood (ml/ml), must be connected externally
    cvCO₂(t) # CO₂ concentration in the deoxygenated blood (ml/ml), must be connected externally


    FDO₂(t) # Fractional concentration of O₂ in the dead space (%)
    FDCO₂(t) # Fractional concentration of CO₂ in the dead space (%)
    FAO₂(t) # Fractional concentration of O₂ in the alveolar space (%)
    FACO₂(t) # Fractional concentration of CO₂ in the alveolar space (%)

    qpp(t) # Pulmonary blood flow (ml/s)
    qps(t) # Pulmonary shunted blood flow (ml/s)
    Vpp(t) # Pulmonary peripheral blood volume (ml)

    cppO₂(t) # O₂ concentration in the pulmonary blood (ml/ml)
    cppCO₂(t) # CO₂ concentration in the pulmonary blood (ml/ml)
    XppO₂(t) # O₂ partial pressure
    XppCO₂(t) # CO₂ partial pressure
    pppO₂(t) # O₂ partial pressure in pulmonary blood (mmHg)
    pppCO₂(t) # CO₂ partial pressure in pulmonary blood (mmHg)

    p_AO₂(t) # O₂ partial pressure in the alveolar space (mmHg)
    p_ACO₂(t) # CO₂ partial pressure in the alveolar space (mmHg)

    caO₂(t) # O₂ concentration in the arterial blood (ml/ml) (output)
    caCO₂(t) # CO₂ concentration in the arterial blood (ml/ml) (output)

    VrO₂p(t) # O₂ uptake (ml/s)
    VrCO₂p(t) # CO₂ uptake (ml/s)
    VrO₂d(t) # O₂ diffusion (ml/s)
    VrCO₂d(t) # CO₂ diffusion (ml/s)
    VrO₂(t) # O₂ uptake (ml/s)
    VrCO₂(t) # CO₂ uptake (ml/s)

    XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
    XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
    paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
    paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

    SaO₂(t) # O₂ saturation in the arterial blood (%)
  end
  @equations begin
    #### Conservation of mass equations
    D(FDO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FIO₂ - FDO₂), Vr_A * (FDO₂ - FAO₂)) / V_D
    D(FDCO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FICO₂ - FDCO₂), Vr_A * (FDCO₂ - FACO₂)) / V_D

    D(FAO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁))

    D(FACO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂))

    Vpp ~ (Vpa - v0pa) * (1 - _sh) + v0pa


    # #### Dissociation equations
    cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
    XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
    cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
    XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

    # #### Instantaneous equilibrium equations
    p_AO₂ ~ pppO₂
    p_ACO₂ ~ pppCO₂

    #### Gas fraction to partial pressure relationships
    p_AO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
    p_ACO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

    #### Mixing between capillary and shunted blood
    caO₂ ~ (1 - _sh) * cppO₂ + _sh * cvO₂
    caCO₂ ~ (1 - _sh) * cppCO₂ + _sh * cvCO₂

    qpp ~ qpa * (1 - _sh)
    qps ~ qpa * _sh

    #### O₂ saturation in arterial blood
    XaO₂ ~ max(((caO₂ / _CₛₐₜO₂) / (1 - (caO₂ / _CₛₐₜO₂))),1e-8)^_h₁
    XaCO₂ ~ max(((caCO₂ / _CₛₐₜCO₂) / (1 - (caCO₂ / _CₛₐₜCO₂))),1e-8)^_h₂

    paO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) + _K₁*XaO₂*_β₂ - _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₂ + 2*_K₂*XaCO₂*_α₂*_β₁)

    paCO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) - _K₁*XaO₂*_β₂ + _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₁ + 2*_K₁*XaO₂*_α₁*_β₂)

    SaO₂ ~ (caO₂ - paO₂ * _sol_O₂) / (_Hgb * _Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
  end
end

"""
Need to decouple the equations:
"""

D(FAO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁))

D(FACO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂))

####Try a substitution

"""
O₂ First
"""

eq = D(FAO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*((-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂)))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*((-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*D(FAO₂)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*D(FAO₂)) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂)))) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁))

solution = solve_for(eq, D(FAO₂))[1]

println(solution)

#### Result

D(FAO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂) + ((-(ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + ((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))) / (-1 + (-_K*(((_K^2)*((-_p_ws + _pₐₜₘ)^2)*_β₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + (-(_K^2)*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))

"""
CO₂ Second
"""

eq = D(FACO₂) ~ (-((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + (-((_CₛₐₜCO₂ * (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂)/(1 + (((FACO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₂ * (FAO₂ * (_pₐₜₘ - _p_ws))) / (_K₂ * (1 + _α₂ * (FAO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₂))) - cvCO₂)*qpp - Vpp*_K*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))*((FACO₂*((-_p_ws + _pₐₜₘ)^2)*_β₂*((-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁)))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂) + (-FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*((-((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + (-((_CₛₐₜO₂ * (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁)/(1 + (((FAO₂ * (_pₐₜₘ - _p_ws)) * (1 + _β₁ * (FACO₂ * (_pₐₜₘ - _p_ws))) / (_K₁ * (1 + _α₁ * (FACO₂ * (_pₐₜₘ - _p_ws))))))^(1/_h₁))) - cvO₂)*qpp + Vpp*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))*((FAO₂*((-_p_ws + _pₐₜₘ)^2)*_β₁*D(FACO₂)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁) + ((-1 - FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*((-_p_ws + _pₐₜₘ)^2)*_α₁*D(FACO₂)) / (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)^2)*_K₁)))*_K) / V_A)) / (-1 + ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*Vpp*_K*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*((1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_β₁)*FAO₂*(-_p_ws + _pₐₜₘ)) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*_K₁))^(1 / _h₁))))) / ((1 + FACO₂*(-_p_ws + _pₐₜₘ)*_α₁)*V_A*_K₁)))) / (((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)^2)*_K₂)))*_K) / V_A)) / (-1 + (-(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*Vpp*(_K^2)*(-_p_ws + _pₐₜₘ)*((-_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))) / (_h₂*((1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))^2)) + (_CₛₐₜCO₂*(((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(-1 + 1 / _h₂))) / (_h₂*(1 + ((FACO₂*(1 + FAO₂*(-_p_ws + _pₐₜₘ)*_β₂)*(-_p_ws + _pₐₜₘ)) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*_K₂))^(1 / _h₂))))) / ((1 + FAO₂*(-_p_ws + _pₐₜₘ)*_α₂)*V_A*_K₂))

solution = solve_for(eq, D(FACO₂))[1]

println(solution)

#### Result

D(FACO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂) - _K*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*(((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-(ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp)) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))) / (-1 + ((_K^2)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*((_K*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-_K*((-_p_ws + _pₐₜₘ)^2)*_β₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))


"""

Take 2

"""

@mtkmodel LungGasExchange begin
  @parameters begin
    _FIO₂ = FIO₂ # Fractional concentration of O₂ in the inspired air (%)
    _FICO₂ = FICO₂ # Fractional concentration of CO₂ in the inspired air (%)
    _sh = sh # Shunt fraction
    _v0pa = v0pa # Unstressed volume of the pulmonary arterial blood (ml)
    _K = K # Proportionality constant that allows convertion of volumes from body temperature pressure saturated to standard temperature pressure dry conditions
    _CₛₐₜO₂ = CₛₐₜO₂ # O₂ saturation constant (ml/ml)
    _CₛₐₜCO₂ = CₛₐₜCO₂ # CO₂ saturation constant (ml/ml)
    _h₁ = h₁ # O₂ saturation exponent
    _h₂ = h₂ # CO₂ saturation exponent
    _K₁ = K₁ # O₂ saturation constant (mmHg)
    _K₂ = K₂ # CO₂ saturation constant (mmHg)
    _α₁ = α₁ # O₂ saturation constant (/mmHg)
    _α₂ = α₂ # CO₂ saturation constant (/mmHg)
    _β₁ = β₁ # O₂ saturation constant (/mmHg)
    _β₂ = β₂ # CO₂ saturation constant (/mmHg)
    _pₐₜₘ = pₐₜₘ # Absolute Atmospheric pressure (mmHg)
    _p_ws = p_ws # Water vapor pressure (mmHg)
    _sol_O₂ = sol_O₂ # Solubility of O₂ in blood (ml O₂/ml blood/mmHg)
    _Hgb_O₂_binding = Hgb_O₂_binding # Hemoglobin O₂ binding constant (ml O₂/ml blood/mmHg)
    _Hgb = Hgb # Hemoglobin concentration (g/dl)
    _DlO₂ = DlO₂ # Diffusion coefficient of O₂ in blood (ml O₂/ml blood/mmHg)
    _DlCO₂ = DlCO₂ # Diffusion coefficient of CO₂ in blood (ml CO₂/ml blood/mmHg)
  end
  @variables begin
    #### Connected from Lung Model
    Vrᵢₙ(t) # Flow rate (ml/s), must be connected externally
    Vr_A(t) # Alveolar flow rate (ml/s), must be connected externally
    V_D(t) # Dead space volume (ml), must be connected externally
    V_A(t) # Alveolar volume (ml), must be connected externally

    #### Connected from Pulmonary Artery
    qpa(t) # Pulmonary blood flow (ml/s), must be connected externally
    Vpa(t) # Pulmonary arterial blood volume (ml), must be connected externally
    cvO₂(t) # O₂ concentration in the deoxygenated blood (ml/ml), must be connected externally
    cvCO₂(t) # CO₂ concentration in the deoxygenated blood (ml/ml), must be connected externally


    FDO₂(t) # Fractional concentration of O₂ in the dead space (%)
    FDCO₂(t) # Fractional concentration of CO₂ in the dead space (%)
    FAO₂(t) # Fractional concentration of O₂ in the alveolar space (%)
    FACO₂(t) # Fractional concentration of CO₂ in the alveolar space (%)

    qpp(t) # Pulmonary blood flow (ml/s)
    qps(t) # Pulmonary shunted blood flow (ml/s)
    Vpp(t) # Pulmonary peripheral blood volume (ml)

    cppO₂(t) # O₂ concentration in the pulmonary blood (ml/ml)
    cppCO₂(t) # CO₂ concentration in the pulmonary blood (ml/ml)
    XppO₂(t) # O₂ partial pressure
    XppCO₂(t) # CO₂ partial pressure
    pppO₂(t) # O₂ partial pressure in pulmonary blood (mmHg)
    pppCO₂(t) # CO₂ partial pressure in pulmonary blood (mmHg)

    p_AO₂(t) # O₂ partial pressure in the alveolar space (mmHg)
    p_ACO₂(t) # CO₂ partial pressure in the alveolar space (mmHg)

    caO₂(t) # O₂ concentration in the arterial blood (ml/ml) (output)
    caCO₂(t) # CO₂ concentration in the arterial blood (ml/ml) (output)

    VrO₂p(t) # O₂ uptake (ml/s)
    VrCO₂p(t) # CO₂ uptake (ml/s)
    VrO₂d(t) # O₂ diffusion (ml/s)
    VrCO₂d(t) # CO₂ diffusion (ml/s)
    VrO₂(t) # O₂ uptake (ml/s)
    VrCO₂(t) # CO₂ uptake (ml/s)

    XaO₂(t) # O₂ saturation in the arterial blood (ml/ml)
    XaCO₂(t) # CO₂ saturation in the arterial blood (ml/ml)
    paO₂(t) # O₂ partial pressure in the arterial blood (mmHg)
    paCO₂(t) # CO₂ partial pressure in the arterial blood (mmHg)

    SaO₂(t) # O₂ saturation in the arterial blood (%)
  end
  @equations begin
    #### Conservation of mass equations
    D(FDO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FIO₂ - FDO₂), Vr_A * (FDO₂ - FAO₂)) / V_D
    D(FDCO₂) ~ ifelse(Vrᵢₙ >= 0, Vrᵢₙ * (_FICO₂ - FDCO₂), Vr_A * (FDCO₂ - FACO₂)) / V_D

    D(FAO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂) + ((-(ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + ((ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂))*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))) / (-1 + (-_K*(((_K^2)*((-_p_ws + _pₐₜₘ)^2)*_β₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A) + (-(_K^2)*(1 - (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*(((-1 - (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))*FAO₂*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)*(-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / ((-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))

    D(FACO₂) ~ (-((-ifelse(Vrᵢₙ >= 0, (-FACO₂ + FDCO₂)*Vr_A, 0) - _K*(qpp*((-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂)) + cvCO₂) - _K*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*(((ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*FACO₂) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-(ifelse(Vrᵢₙ >= 0, (-FAO₂ + FDO₂)*Vr_A, 0) + _K*qpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁)) + cvO₂))*((-_p_ws + _pₐₜₘ)^2)*_β₂*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp)) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))) / (-1 + ((_K^2)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*((_K*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*((1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)^2)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A) + (-_K*((-_p_ws + _pₐₜₘ)^2)*_β₂*Vpp*FACO₂*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))*(((-1 + (_p_ws - _pₐₜₘ)*_β₁*FACO₂)*((-_p_ws + _pₐₜₘ)^2)*_α₁*FAO₂) / (_K₁*((1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)^2)) + (((-_p_ws + _pₐₜₘ)^2)*_β₁*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*(-1 + (_K*(-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*Vpp*((-_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))) / (_h₁*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))^2)) + (_CₛₐₜO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(-1 + 1 / _h₁))) / (_h₁*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₁*FACO₂)*FAO₂) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)))^(1 / _h₁))))) / (_K₁*(1 + (-_p_ws + _pₐₜₘ)*_α₁*FACO₂)*V_A))*V_A))*Vpp) / ((-1 + ((_K^2)*(-_p_ws + _pₐₜₘ)*(-1 + (_p_ws - _pₐₜₘ)*_β₂*FAO₂)*((_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))) / (_h₂*(1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) + (-_CₛₐₜCO₂*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(-1 + 1 / _h₂))*((((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))) / (_h₂*((1 + (((-_p_ws + _pₐₜₘ)*(1 + (-_p_ws + _pₐₜₘ)*_β₂*FAO₂)*FACO₂) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)))^(1 / _h₂))^2)))*Vpp) / (_K₂*(1 + (-_p_ws + _pₐₜₘ)*_α₂*FAO₂)*V_A))*V_A))

    Vpp ~ (Vpa - v0pa) * (1 - _sh) + v0pa


    # #### Dissociation equations
    cppO₂ ~ _CₛₐₜO₂ * (XppO₂)^(1/_h₁)/(1 + (XppO₂)^(1/_h₁))
    XppO₂ ~ pppO₂ * (1 + _β₁ * pppCO₂) / (_K₁ * (1 + _α₁ * pppCO₂))
    cppCO₂ ~ _CₛₐₜCO₂ * (XppCO₂)^(1/_h₂)/(1 + (XppCO₂)^(1/_h₂))
    XppCO₂ ~ pppCO₂ * (1 + _β₂ * pppO₂) / (_K₂ * (1 + _α₂ * pppO₂))

    # #### Instantaneous equilibrium equations
    p_AO₂ ~ pppO₂
    p_ACO₂ ~ pppCO₂

    #### Gas fraction to partial pressure relationships
    p_AO₂ ~ FAO₂ * (_pₐₜₘ - _p_ws)
    p_ACO₂ ~ FACO₂ * (_pₐₜₘ - _p_ws)

    #### Mixing between capillary and shunted blood
    caO₂ ~ (1 - _sh) * cppO₂ + _sh * cvO₂
    caCO₂ ~ (1 - _sh) * cppCO₂ + _sh * cvCO₂

    qpp ~ qpa * (1 - _sh)
    qps ~ qpa * _sh

    #### O₂ saturation in arterial blood
    XaO₂ ~ max(((caO₂ / _CₛₐₜO₂) / (1 - (caO₂ / _CₛₐₜO₂))),1e-8)^_h₁
    XaCO₂ ~ max(((caCO₂ / _CₛₐₜCO₂) / (1 - (caCO₂ / _CₛₐₜCO₂))),1e-8)^_h₂

    paO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) + _K₁*XaO₂*_β₂ - _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₂ + 2*_K₂*XaCO₂*_α₂*_β₁)

    paCO₂ ~ (-1 + √(1 + 2*_K₁*XaO₂*_β₂ + 2*_K₂*XaCO₂*_β₁ + (_K₁^2)*(XaO₂^2)*(_β₂^2) - 2*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₁*_β₂ + 4*_K₁*_K₂*XaCO₂*XaO₂*_α₂*_β₁ - 2*_K₁*_K₂*XaCO₂*XaO₂*_β₁*_β₂ + (_K₂^2)*(XaCO₂^2)*(_β₁^2) + 2(_K₁^2)*_K₂*XaCO₂*(XaO₂^2)*_α₁*_α₂*_β₂ + 2*_K₁*(_K₂^2)*(XaCO₂^2)*XaO₂*_α₁*_α₂*_β₁ + (_K₁^2)*(_K₂^2)*(XaCO₂^2)*(XaO₂^2)*(_α₁^2)*(_α₂^2)) - _K₁*XaO₂*_β₂ + _K₂*XaCO₂*_β₁ + _K₁*_K₂*XaCO₂*XaO₂*_α₁*_α₂) / (2*_β₁ + 2*_K₁*XaO₂*_α₁*_β₂)

    SaO₂ ~ (caO₂ - paO₂ * _sol_O₂) / (_Hgb * _Hgb_O₂_binding) * 100 # O₂ saturation in arterial blood (%)
  end
end

"""
Take 3
"""

eq = D(FAO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDO₂ - FAO₂),0) - _K * (qpp * (cppO₂ - cvO₂) - Vpp * D(cppO₂))) / V_A

sol2 = solve_for(eq, D(cppO₂))[1]

(((-ifelse(Vrᵢₙ(t) >= 0, (-FAO₂(t) + FDO₂(t))*Vr_A(t), 0) + _K*qpp(t)*(cppO₂(t) - cvO₂(t))) / V_A(t) + Differential(t)(FAO₂(t)))*V_A(t)) / (_K*Vpp(t))

D(FACO₂) ~ (ifelse(Vrᵢₙ >= 0, Vr_A * (FDCO₂ - FACO₂),0) - _K * (qpp * (cppCO₂ - cvCO₂) + _K * Vpp * D(cppCO₂))) / V_A