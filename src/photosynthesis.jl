abstract type AbstractPhotosynthesisModel end

@with_kw struct Model_Photosynthesis{FT} <: AbstractPhotosynthesisModel{FT}
  "initial slope of the light response curve to assimilation rate, (i.e., quantum efficiency; `μmol CO2 [μmol PAR]⁻¹`)`"
  α::FT = Param(0.06, bounds=(0.01, 0.10))

  "initial slope of the CO2 response curve to assimilation rate, (i.e., carboxylation efficiency; `μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹`)"
  η::FT = Param(0.04, bounds=(0.01, 0.07))

  "carbon saturated rate of photosynthesis at 25 °C, `μmol m⁻² s⁻¹`"
  VCmax25::FT = Param(50.00, bounds=(5.00, 120.00))

  "photoperiod constraint"
  d_pc::FT = Param(2.0, bounds=(0.0, 5.0))
end


"""
    photosynthesis(Tavg::T, Rs::T, VPD::T, LAI::T, Ca=380.0; par)

# Example
```julia
# GPP, Gc_w = photosynthesis(Tavg, Rs, VPD, LAI, Ca; par)
```
"""
function photosynthesis(Tavg::T, Rs::T, VPD::T, LAI::T,
  Pa=atm, Ca=380.0, PC=1.0; par::Param_PMLV2) where {T<:Real}
  (; kQ, VCmax25, α, η, g1, d_pc, D0, VPDmin, VPDmax) = par

  PAR = 0.45 * Rs # W m-2, taken as 0.45 time of solar radiation
  PAR_mol = PAR * 4.57 # 1 W m-2 = 4.57 umol m-2 s-1

  Vm = VCmax25 * T_adjust_Vm25(Tavg) * PC^d_pc # * data$dhour_norm^2 
  Am = Vm # 认为最大光合速率 = 最大羧化能力

  P1 = Am * α * η * PAR_mol
  P2 = Am * α * PAR_mol
  P3 = Am * η * Ca
  P4 = α * η * PAR_mol * Ca

  ## canopy conductance in (mol m-2 s-1)
  Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
    kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4))) # [umol m-2 s-1]
  Ag = Ags  # gross assimilation rate in [umol m-2 s-1]
  Ag = Ag * f_VPD_Zhang2019(VPD, VPDmin, VPDmax)

  GPP = Ag * 86400 / 10^6 * 12 # [umol m-2 s-1] to [g C m-2 d-1]

  f_VPD_gc = 1.0 / (1.0 + VPD / D0) # Leuning f_vpd
  Gc = g1 * Ag / Ca * f_VPD_gc # canopy conductance for carbon, [umol m-2 s-1] / [umol mol-1] = [mol m-2 s-1]

  ## Convert from [mol m-2 s-1] to m s-1
  Gc = Gc * 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3)) # Gc = Gc * mol2m(Tavg, Pa)
  Gc = max(Gc, 1e-6)

  Gc_w = Gc * 1.6 # g_water  = 1.6 * g_CO2 (mol m-2 s-1), canopy conductance for water
  GPP, Gc_w
end


# 最大羧化能力温度调节函数
# V_m = Vm_25 * T_adjust_Vm25
# `T` in -100:100, `T_adjust_Vm25` always < 0.92
function T_adjust_Vm25(Tavg::T)::T where {T<:Real}
  a = 0.031
  b = 0.115
  exp(a * (Tavg - 25.0)) / (1.0 + exp(b * (Tavg - 41.0))) # Gan2018, Eq. A5
end



"""
    mol2m(x, Tavg, Pa=atm)
    mol2m_rong2018(x, Tavg, Pa=atm)

Convert from mol m-2 s-1 to m s-1, g = g_m * mol2m(Tavg)

## Arguments
- `Tavg`: degC
- `Pa`: kPa

# Reference
1. Monteith, 2013, Principles of Environmental Physics, Eq. 3.14
1. CLM5 TechNote, Eq. 9.1
"""
mol2m(x, Tavg, Pa=atm) = x * R * (Tavg + K0) / (Pa * 1000)

mol2m_rong2018(x, Tavg, Pa=atm) = x * 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3))
