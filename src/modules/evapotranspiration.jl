"""
    Evapotranspiration_PML{FT<:Real} <: AbstractEvapotranspirationModel{FT}

# Fields
$(TYPEDFIELDS)
"""
@bounds @with_kw mutable struct Evapotranspiration_PML{FT<:Real} <: AbstractEvapotranspirationModel{FT}
  "extinction coefficients for available energy"
  kA::FT = 0.70 | (0.50, 0.9)

  "Specific leaf storage, van Dijk, A.I.J.M, 2001, Eq2"
  S_sls::FT = 0.1 | (0.01, 1.0)

  "Canopy cover fraction related parameter"
  fER0::FT = 0.1 | (0.01, 0.5)

  "canopy height, `[m]`"
  hc::FT = 1.0 | (0.01, 20.0)

  # LAIref::FT = 4.0   | (1.0, 6.0)      # 
  # frame::Integer = 10.0  | (6.0, 14.0) # 8-day moving window
end


"""
    PMLV2 (Penman–Monteith–Leuning Version 2) Evapotranspiration model

# Arguments

- `Prcp` : mm/d
- `Tavg` : degC
- `Rs`   : W m-2
- `Rn`   : W m-2
- `VPD`  : W m-2
- `U2`   : m/s
- `LAI`  : m2 m-2
- `Pa`   : kPa
- `Ca`   : ppm, default `380`
- `Ω`    : clamping index, default is `1.0`

# Examples
```julia
```

# References
1. Gan Rong, 2018, Ecohydrology
2. Zhang Yongqiang, 2019, RSE
3. Kong Dongdong, 2019, ISPRS
"""
function evapotranspiration(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel,
  stomatal::AbstractStomatalModel,
  Prcp::T, Tavg::T, Rs::T, Rn::T, VPD::T, U2::T, LAI::T,
  Pa=atm, Ca=380.0, PC=1.0,
  Ω::T=T(1.0);
  r::Union{Nothing,interm_PML}=nothing) where {T<:Real}

  (; fER0, S_sls) = evap

  isnothing(r) && (r = interm_PML{T}())
  
  ## Intercepted Evaporation (Ei)
  λ::T = cal_lambda(Tair) # [MJ kg-1]
  γ::T = cal_gamma(Tair)
  Δ::T = cal_slope(Tair)

  Ei = cal_Ei_Dijk2021(Prcp, LAI, fER0, S_sls)
  Pi = Prcp - Ei # 应扣除这一部分消耗的能量

  Rn = Rn - MJ2W(λ * Pi)

  τ = exp(-kA * LAI)
  Rns = (1 - τ) * Rn
  Rnc = τ * Rn

  rs = leaf_conductance(model, met) # coupling GPP model
  Ec, Ecr, Eca, ra = ET0_Monteith65(Rnc, Tair, VPD, U2, Pa;
    rs, hc, z_wind=2.0)

  # 水面蒸发
  Eeq_water = Δ / (Δ + γ) * Rn |> x -> W2mm(x, λ)
  Evp_water = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ # [mm]
  ET_water = Eeq_water + Evp_water

  ## TODO: 补充冰面蒸发的计算
  Es_eq = Δ / (Δ + γ) * Rns |> x -> W2mm(x, λ) # 土壤均衡蒸发
  # r.Es_eq = r.Eeq * exp(-0.9 * LAI) # Ka = 0.9, insensitive param
  # r
  (GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, ET_water, rs, ra)
end


function leaf_conductance(
  photo::AbstractPhotosynthesisModel,
  stomatal::AbstractStomatalModel,
  Tavg, Rs, VPD, LAI, Pa, Ca, PC;)

  Ag, Rd = photosynthesis(photo, Tavg, Rs, VPD, LAI, Pa, Ca, PC;)
  GPP = umol2gC(Ag)

  gs = stomatal_conductance(stomatal, Ag, Rd, VPD, Ca) # [mol m-2 s-1], 0.1~0.4
  rs = 1 / mol2m(gs, Tavg, Pa)
  rs
end
