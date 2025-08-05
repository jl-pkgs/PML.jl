export LandModel, OverUnderCanopy
export AbstractWaterConsGPPModel, AbstractStomatalModel, AbstractPhotosynthesisModel, AbstractEvapotranspirationModel
export AirLayer, update!


using Parameters
import FieldMetadata: @bounds, bounds, @units, units

# abstract type AbstractModel{FT} end

# 水分限制-GPP
abstract type AbstractWaterConsGPPModel{FT} <: AbstractModel{FT} end

# 气孔导度
abstract type AbstractStomatalModel{FT} <: AbstractModel{FT} end

# 光合
abstract type AbstractPhotosynthesisModel{FT} <: AbstractModel{FT} end

# 蒸发
abstract type AbstractEvapotranspirationModel{FT} <: AbstractModel{FT} end



@units @with_kw mutable struct AirLayer{FT<:Real}
  Prcp::FT = 0.0 | "mm d⁻¹"
  Tavg::FT = 20.0 | "℃"
  Rs::FT = 200.0 | "W m⁻²"
  Rn::FT = 150 | "W m⁻²"
  VPD::FT = 2.0 | "kPa"
  U2::FT = 1.0 | "m s⁻¹"
  Pa::FT = atm | "kPa"
  Ca::FT = 380.0 | "μmol mol⁻¹" # air CO2 content
  PC::FT = 1.0 | "-" # photoperiod constaint

  # auxiliary vairables
  λ::FT = cal_lambda(Tavg) | "MJ kg⁻¹"
  γ::FT = cal_gamma(Tavg, Pa) | "kPa K⁻¹"
  Δ::FT = cal_slope(Tavg) | "kPa K⁻¹"
  β::FT = γ / Δ | "-" # Bowen ratio, H / LE, γ / Δ
  ρₐ::FT = cal_rho_a(Tavg, Pa) | "kg m⁻³"
end


function update!(air::AirLayer{FT}; Tavg::FT, Pa::FT, kw...) where {FT<:Real}
  λ::FT = cal_lambda(Tavg) # "MJ kg⁻¹"
  γ::FT = cal_gamma(Tavg, Pa) # "kPa K⁻¹"
  Δ::FT = cal_slope(Tavg) # "kPa K⁻¹"
  β::FT = γ / Δ # "-" # Bowen ratio, H / LE, γ / Δ
  ρₐ::FT = cal_rho_a(Tavg, Pa) # "kg m⁻³"
  @pack! air = Tavg, Pa, λ, γ, Δ, β, ρₐ, kw...
end




@with_kw mutable struct LandModel{FT} <: AbstractModel{FT}
  # atmospheric_forcing::AtmosphericForcing
  evap::AbstractEvapotranspirationModel{FT}
  photo::AbstractPhotosynthesisModel{FT}
  stomatal::AbstractStomatalModel{FT}
end
