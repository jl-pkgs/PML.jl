using Dates


"""
任一点，任一时刻的太阳高度角
"""
function angle_SunElevation(lat::Real, time_local::DateTime)
  ψ = deg2rad(lat)  # 纬度转弧度
  J = dayofyear(time_local)
  σ = angle_SolarDeclination(J)
  dh = hour(time_local) + minute(time_local) / 60 - 12.0
  ω = deg2rad(dh * 15) # 时角，上午为负，下午为正

  sinH = cos(ψ) * cos(σ) * cos(ω) + sin(ψ) * sin(σ)
  return asin(sinH)
end


"Solar Declination Angle (σ), 黄赤交角(太阳赤纬角)"
angle_SolarDeclination(J::Integer) = 0.409 * sin(2pi / 365 * J - 1.39)  # Allen 1998, Eq. 24


"""
    叶片角度分布--球体偏离度

球体(0)，水平叶片(1)，垂直叶片(-1)

- F1: 倾角 0-30°叶片比例
- F2: 倾角30-60°叶片比例
- F3: 倾角60-90°叶片比例
"""
function cal_χₗ(F1::T, F2::T, F3::T) where {T<:Real}
  χₗ = sign(0.5 - F3) / 2 * (abs(0.134 - F1) + abs(0.366 - F2) + abs(0.5 - F3)) # Bonan 2019, Eq. 2.17
  if χₗ > 0.6 || χₗ < -0.4
    @warn "χₗ is out [-0.4, 0.6]"
    χₗ = clamp(χₗ, T(-0.4), T(0.6))
  end
  return χₗ
end


"""
- `χₗ`: 球体偏离度。球体(0)，水平叶片(1)，垂直叶片(-1)。
- `Z`: 太阳天顶角，[0, pi/2]，中午为0。sinh = cosZ
"""
function cal_G(χₗ::T, cosZ::T) where {T}
  Φ1 = 0.5 - 0.633 * χₗ - 0.330 * χₗ^2
  Φ2 = 0.877 * (1 - 2 * Φ1)
  # This equation is restricted to 0:4 < χℓ < 0:6
  Φ1 + Φ2 * cosZ # Bonan 2019, Eq. 14.31
end
