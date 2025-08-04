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


@with_kw mutable struct LandModel{FT} <: AbstractModel{FT}
  # atmospheric_forcing::AtmosphericForcing
  stomatal::AbstractStomatalModel{FT}
  photosynthesis::AbstractPhotosynthesisModel{FT}
  evapotranspiration::AbstractEvapotranspirationModel{FT}
  watercons_GPP::AbstractWaterConsGPPModel{FT}
end
