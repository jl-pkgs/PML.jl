# export struct2vec, struct2tuple
export round2;
export movmean2, nanmean2
export map_df_tuple

# rounded_data = NamedTuple((field => round(value) for (field, value) in data))
round2(x::NamedTuple, digits=3; kw...) = map(val -> round(val; digits), x)


function nanmean2(x::T1, y::T2) where {T1<:Real, T2<:Real}
  T = promote_type(T1, T2)
  if isnan(x)
    T(y)
  elseif isnan(y)
    T(x)
  else
    (x + y) / 2
  end
end


# function struct2vec(x)
#   keys = fieldnames(typeof(x))
#   vals = [getfield(x, key) for key in keys]
#   vals
# end

# function struct2tuple(x)
#   keys = fieldnames(typeof(x))
#   vals = [getfield(x, key) for key in keys]
#   (; zip(keys, vals)...)
# end
weighted_mean(x::AbstractVector, w::AbstractVector) = sum(x .* w) / sum(w)


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
