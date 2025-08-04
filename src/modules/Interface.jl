export LandModel
export AbstractWaterConsGPPModel, AbstractStomatalModel, AbstractPhotosynthesisModel, AbstractEvapotranspirationModel


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
  Prcp::FT = 0.0 | "mm"
  Tavg::FT = 20.0 | "℃"
  Rs::FT = 200.0 | "W m⁻²"
  Rn::FT = 150 | "W m⁻²"
  VPD::FT = 2.0 | "kPa"
  U2::FT = 1.0 | "m s⁻¹"
  Pa::FT = atm | "kPa"
  Ca::FT = 380.0 | "μmol mol⁻¹" # air CO2 content
  PC::FT = 1.0 | "-" # photoperiod constaint
end


@units @with_kw mutable struct CanopyLayer{FT<:Real}
  LAI::FT = 2.0 | "m² m⁻²" # leaf area index
  Ω::FT = 1.0 | ""  # clamping index, default is `1.0`
end



@with_kw mutable struct LandModel{FT} <: AbstractModel{FT}
  # atmospheric_forcing::AtmosphericForcing
  stomatal::AbstractStomatalModel{FT}
  photosynthesis::AbstractPhotosynthesisModel{FT}
  # evapotranspiration::AbstractEvapotranspirationModel{FT}
end
