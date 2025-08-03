using Parameters, DataFrames
import FieldMetadata: @bounds, bounds


abstract type AbstractModel{FT} end

# 水分限制-GPP
abstract type AbstractWaterConsGPPModel{FT} <: AbstractModel{FT} end

# 气孔导度
abstract type AbstractStomatalModel{FT} <: AbstractModel{FT} end

# 光合
abstract type AbstractPhotosynthesisModel{FT} <: AbstractModel{FT} end

# 蒸发
abstract type AbstractEvapotranspirationModel{FT} <: AbstractModel{FT} end

struct_values(s::T) where {T} = [getfield(s, name) for name in fieldnames(T)]
struct_names(::T) where {T} = collect(fieldnames(T))


@with_kw struct LandModel{FT}
  # atmospheric_forcing::AtmosphericForcing
  stomatal::AbstractStomatalModel{FT}
  photosynthesis::AbstractPhotosynthesisModel{FT}
  evapotranspiration::AbstractEvapotranspirationModel{FT}
  watercons_GPP::AbstractWaterConsGPPModel{FT}
end


function get_params(mod::AbstractModel; kw...)
  DataFrame(; kw..., name=struct_names(mod), value=struct_values(mod), bound=collect(bounds(mod)))
end

function get_params(model::LandModel)
  modules = fieldnames(typeof(model))
  params = map(f -> get_params(getfield(model, f); mod=f), modules)
  vcat(params...)
end


## 设置模型参数
function update!(model::LandModel, key::Symbol, value::T;
  pars_proxy::Union{Nothing,DataFrame}=nothing) where {T<:Real}
  isnothing(pars_proxy) && (pars_proxy = get_params(model))

  mod_name = pars_proxy |> d -> d[d.name.==key, :mod]
  @assert length(mod_name) == 1 "Duplicated parameter name found!"

  mod = getfield(model, mod_name[1])
  setfield!(mod, key, value)
end


"""
    update!(model::LandModel, values::Vector{T}, keys::Vector{Symbol})

# Examples
```julia
update!(model, keys, values)
```
模型参数名不能有重复！
"""
function update!(model::LandModel, keys::Vector{Symbol}, values::Vector{T};
  pars_proxy::Union{Nothing,DataFrame}=nothing) where {T<:Real}
  isnothing(pars_proxy) && (pars_proxy = get_params(model))

  for (key, value) in zip(keys, values)
    update!(model, key, value)
  end
end
