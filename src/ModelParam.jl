export ModelParams, AbstractModel, update!


import FieldMetadata: @bounds, bounds, @units, units

abstract type AbstractModel{FT} end

function unlist(list::Vector)
  res = []
  for x in list
    isa(x, Vector) ? append!(res, unlist(x)) : push!(res, x)
  end
  return map(x -> x, res)
end


_fieldname(x, field::Symbol) = field

"递归获取模型参数属性"
function get_ModelParamRecursive(x::T, path=Vector{Symbol}(); fun=bounds) where {FT,T<:AbstractModel{FT}}
  fileds = fieldnames(T) |> collect
  map(field -> begin
      value = getfield(x, field)
      _path = [path..., field]
      r = fun(x, field)
      typeof(value) == FT ? (_path => r) : get_ModelParamRecursive(value, _path; fun)
    end, fileds) |> unlist
end


function ModelParams(model::AbstractModel)
  _names = get_ModelParamRecursive(model; fun=_fieldname)
  _bounds = get_ModelParamRecursive(model; fun=bounds)
  _values = get_ModelParamRecursive(model; fun=getfield)
  _units = get_ModelParamRecursive(model; fun=units)

  # 返回NamedTuple向量，符合Tables接口
  [(path=first(n), name=last(n), value=last(v),
    unit=last(u), bound=last(b))
   for (n, v, u, b) in zip(_names, _values, _units, _bounds)]
end


function update!(model::AbstractModel{FT}, names::Vector{Symbol}, values::Vector{FT},
  ; params::Union{Nothing,Vector{<:NamedTuple}}=nothing) where {FT}
  isnothing(params) && (params = ModelParams(model))

  for (name, value) in zip(names, values)
    rows = filter(row -> row.name == name, params)
    @assert length(rows) == 1 "Duplicated parameters are not allowed!"
    update!(model, rows[1].path, value)
  end
end

function update!(model::AbstractModel{FT}, path::Vector{Symbol}, value::FT) where {FT}
  if length(path) == 1
    setfield!(model, path[1], value)
  elseif length(path) > 1
    submodel = getfield(model, path[1])
    update!(submodel, path[2:end], value)
  end
end
