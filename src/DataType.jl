@with_kw mutable struct interm_PML{T}
  ET::T = 0.0
  GPP::T = 0.0
  Ec::T = 0.0
  Ecr::T = 0.0
  Eca::T = 0.0

  Ei::T = 0.0
  Pi::T = 0.0
  Es_eq::T = 0.0
  Eeq::T = 0.0
  ET_water::T = 0.0

  Ga::T = 0.0
  Gc_w::T = 0.0

  fval_soil::T = 0.0
  Es::T = 0.0
end


@with_kw mutable struct output_PML{T}
  n::Integer
  ET::Vector{T} = zeros(T, n)
  GPP::Vector{T} = zeros(T, n)
  Ec::Vector{T} = zeros(T, n)
  Ecr::Vector{T} = zeros(T, n)
  Eca::Vector{T} = zeros(T, n)

  Ei::Vector{T} = zeros(T, n)
  Pi::Vector{T} = zeros(T, n)
  Es_eq::Vector{T} = zeros(T, n)
  Eeq::Vector{T} = zeros(T, n)
  ET_water::Vector{T} = zeros(T, n)

  Ga::Vector{T} = zeros(T, n)
  Gc_w::Vector{T} = zeros(T, n)

  fval_soil::Vector{T} = zeros(T, n)
  Es::Vector{T} = zeros(T, n)
end
# output_PML(;n::Integer) = output_PML{Float64}(;n=n)


Base.getindex(x::interm_PML, key::Union{Symbol,String}) = getfield(x, Symbol(key))

function Base.setindex!(res::output_PML, r::Union{NamedTuple,interm_PML}, t, fields)
  # fields = fieldnames(output_PML)[2:end]
  for i = eachindex(fields)
    field = fields[i]
    x = getfield(res, field)
    x[t] = r[field]
  end
end

# function Base.setindex!(res::output_PML, r::NTuple, t, fields)
#   for i = eachindex(fields)
#     field = fields[i]
#     x = getfield(res, field)
#     x[t] = r[i]
#   end
# end


## DATATYPE CONVERSION ---------------------------------------------------------
function to_mat(res::output_PML{T}) where {T<:Real}
  names = fieldnames(output_PML)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data, names
end

function to_df(res::output_PML{T}) where {T<:Real}
  data, names = to_mat(res)
  DataFrame(data, names)
end

function map_df_tuple(fun::Function, lst::GroupedDataFrame{DataFrame}, args...; kw...)
  n = length(lst)
  _keys = keys(lst)
  map(i -> begin
      d = lst[i]
      key = NamedTuple(_keys[i])
      r = fun(d, args...; kw...)
      (; key..., r...)
    end, 1:n)
end


export interm_PML, output_PML
export to_mat;
export map_df_tuple
