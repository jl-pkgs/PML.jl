using Statistics

# rounded_data = NamedTuple((field => round(value) for (field, value) in data))
round2(x::NamedTuple, digits=3; kw...) = map(val -> round(val; digits), x)


function movmean2(y::AbstractVector{T}, win_left::Integer, win_right::Integer=win_right) where {T<:Real}
  n = length(y)
  z = zeros(Float64, n)

  @inbounds for i in 1:n
    i_beg = max(i - win_left, 1)
    i_end = min(i + win_right, n)

    count = 0    # number
    ∑ = T(0.0)   # sum of values in window
    for j in i_beg:i_end
      isnan(y[j]) && continue # skip missing values
      count += 1
      ∑ += y[j]
    end
    z[i] = count > 0 ? ∑ / count : T(NaN)
  end
  return z
end

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

function replace_miss(x::AbstractVector, miss=NaN)
  x[ismissing.(x)] .= miss
  x
end

function getDataType(x)
  type = eltype(x)
  typeof(type) == Union ? type.b : type
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

# export struct2vec, struct2tuple
export round2;
export movmean2, nanmean2, getDataType, replace_miss
