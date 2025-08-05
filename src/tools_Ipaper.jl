# 用于计算β_Es
function movmean2!(z::AbstractVector{T}, y::AbstractVector{T}, win_left::Integer, win_right::Integer=0) where {T<:Real}
  n = length(y)
  # z = zeros(Float64, n)
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


function movmean2(y::AbstractVector{T}, win_left::Integer, win_right::Integer=0) where {T<:Real}
  z = zeros(Float64, y)
  movmean2!(z, y, win_left, win_right)
end
