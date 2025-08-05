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
function evapotranspiration!(
  output::SpacOutput{T},
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

  ## 需要根据前期 Pi/Es_eq 去计算土壤水限制β，因此无法立即给出Es
  Es_eq = Δ / (Δ + γ) * Rn_s |> x -> W2mm(x, λ) # 土壤均衡蒸发
  @pack! output = GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, ET_water, rs, ra
  return output
end

function evapotranspiration(
  air::AirLayer{T}, canopy::BigLeaf{T},
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T<:Real}

  output = SpacOutput{T}()
  evapotranspiration!(output, air, canopy, evap, photo, stomatal)
end


"""
    PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca; par=param0, frame=3)

# Notes
一个站点的计算。注意，不同植被类型，参数不同。

# Arguments
- `frame`: in 8-days
"""
function evapotranspiration(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  Lai::V,
  Prcp::V, Tavg::V, Rs::V, Rn::V,
  VPD::V, U2::V,
  Pa::V, Ca::V, PC::Union{T,V}=T(1.0);
  frame::Int=3,
  res::Union{Nothing,SpacOutputs}=nothing) where {T<:Real,V<:AbstractVector{T}}

  # 是否开启光周期
  (; use_PC) = photo
  !isa(PC, Vector) && (use_PC = false)

  canopy = BigLeaf{T}() # 先测试大叶模型
  air = AirLayer{T}()
  r = SpacOutput{T}()


  n = length(Prcp)
  res === nothing && (res = SpacOutputs{T}(; n))

  @inbounds for t = 1:n
    _PC = use_PC ? PC[t] : T(1.0)
    update!(air, Prcp[i], Tavg[i], Rs[i], Rn[i], VPD[i], U2[i], Pa[i]; Ca=Ca, PC=_PC)
    canopy.Lai = Lai[i]

    evapotranspiration!(r, air, canopy, evap, photo, stomatal)
    res[i] = r
  end

  res.β_Es = movmean2(res.Pi, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.β_Es, 0, 1)

  res.Es .= res.β_Es .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end
