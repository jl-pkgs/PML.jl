export Evapotranspiration_PML, evapotranspiration


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

# References
1. Gan Rong, 2018, Ecohydrology
2. Zhang Yongqiang, 2019, RSE
3. Kong Dongdong, 2019, ISPRS
"""
function evapotranspiration(
  air::AirLayer{T},
  canopy::BigLeaf{T},
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T<:Real}

  (; Prcp, Tavg, Rn, VPD, U2, Pa) = air
  (; Lai) = canopy
  (; fER0, S_sls, kA, hc) = evap

  ## 水面蒸发
  λ::T = cal_lambda(Tavg) # [MJ kg-1]
  γ::T = cal_gamma(Tavg)
  Δ::T = cal_slope(Tavg)

  Eeq_water = Δ / (Δ + γ) * Rn |> x -> W2mm(x, λ)
  Evp_water = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ # [mm]
  ET_water = Eeq_water + Evp_water

  ## Intercepted Evaporation (Ei)
  Ei = cal_Ei_Dijk2021(Prcp, Lai, fER0, S_sls)
  Pi = Prcp - Ei # 应扣除这一部分消耗的能量
  Rn = Rn - MJ2W(λ * Pi)

  radiative_transfer!(canopy, Rn; kA)
  (; Rn_c, Rn_s) = canopy

  GPP, rs = leaf_conductance(air, canopy, photo, stomatal) # coupling GPP model
  Ec, Ecr, Eca, ra = ET0_Monteith65(Rn_c, Tavg, VPD, U2, Pa;
    rs, hc, z_wind=2.0) # fix this part

  # TODO: 补充冰面蒸发的计算
  Es_eq = Δ / (Δ + γ) * Rn_s |> x -> W2mm(x, λ) # 土壤均衡蒸发
  (; GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, ET_water, rs, ra)
end
# r.Es_eq = r.Eeq * exp(-0.9 * LAI) # Ka = 0.9, insensitive param

function transpiration()
end

function evaporation()
end


function leaf_conductance(
  air::AirLayer{T},
  canopy::AbstractLeaf{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai) = canopy
  Lai <= 0.01 && return T(0.0), T(0.0) # GPP, rs # 无冠层无rs

  Ag, Rd = photosynthesis(photo, Tavg, Rs, VPD, Lai, Ca, PC)
  GPP = umol2gC(Ag)

  gs = stomatal_conductance(stomatal, Ag, Rd, VPD, Ca) # [mol m-2 s-1], 0.1~0.4
  rs = 1 / mol2m(gs, Tavg, Pa)
  GPP, rs
end
