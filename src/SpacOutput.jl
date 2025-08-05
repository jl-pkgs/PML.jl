export SpacOutput, SpacOutputs
export struct2NT

import Base: getindex

function struct2NT(x)
  names = fieldnames(typeof(x))
  NamedTuple{names}(getfield(x, name) for name in names)
end

@with_kw mutable struct SpacOutput{T}
  ET::T = 0.0
  GPP::T = 0.0
  Ec::T = 0.0
  Ecr::T = 0.0
  Eca::T = 0.0

  Ei::T = 0.0 # Interception
  Pi::T = 0.0
  Es_eq::T = 0.0
  # Eeq::T = 0.0 # not used
  ET_water::T = 0.0

  ra::T = 0.0
  rs::T = 0.0 # for water

  β_Es::T = 0.0  # 土壤蒸发限制因子, sum(Pi)/sum(Es_eq)
  β_GPP::T = 0.0 # 植被光合限制因子
  Es::T = 0.0
end

function Base.getindex(x::SpacOutput, names::Vector{Symbol})
  [getfield(x, name) for name in names]
end

function Base.getindex(x::SpacOutput, index::Vector{Int})
  names = fieldnames(typeof(x))[index]
  [getfield(x, name) for name in names]
end



@with_kw mutable struct SpacOutputs{FT}
  ## Original Output
  ntime::Int = 10

  ET::Vector{FT} = zeros(ntime)
  GPP::Vector{FT} = zeros(ntime)
  Ec::Vector{FT} = zeros(ntime)
  Ecr::Vector{FT} = zeros(ntime)
  Eca::Vector{FT} = zeros(ntime)

  Ei::Vector{FT} = zeros(ntime)    # Interception
  Pi::Vector{FT} = zeros(ntime)
  Es_eq::Vector{FT} = zeros(ntime) # forced by Rn_s
  ET_water::Vector{FT} = zeros(ntime)
  # Eeq::Vector{FT} = zeros(ntime) # not used

  ra::Vector{FT} = zeros(ntime) # aerodynamic conductance
  rs::Vector{FT} = zeros(ntime) # for water

  β_Es::Vector{FT} = zeros(ntime)  # 土壤水限制因子, sum(Pi)/sum(Es_eq)
  β_GPP::Vector{FT} = zeros(ntime) # 限制因子, sum(Pi)/sum(Es_eq)
  Es::Vector{FT} = zeros(ntime)
end


function Base.setindex!(res::SpacOutputs, r::SpacOutput, i::Int64)
  fields = fieldnames(SpacOutput)
  @inbounds for f in fields
    getfield(res, f)[i] = getfield(r, f)
  end
  return res
end


## DATATYPE CONVERSION ---------------------------------------------------------
import DataFrames: DataFrame

function to_mat(res::SpacOutputs{T}) where {T<:Real}
  TYPE = typeof(res)
  names = fieldnames(TYPE)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data, names
end

function DataFrame(res::SpacOutputs{T}) where {T<:Real}
  data, names = to_mat(res)
  DataFrame(data, names)
end

export to_mat, DataFrame;
