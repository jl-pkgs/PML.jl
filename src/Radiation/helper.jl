export cal_fsun, cal_τb_sum, cal_τb, cal_τd

using Dates

"""
blackbody radiation

# Arguments
- `Ta`: air temperature, in degC
- `ϵ`: emissivity, default is 1

# Return
- longwave radiation, in `W m-2`
"""
function blackbody(ϵ::T, Ta::T) where {T<:Real}
  σ = 5.67e-8
  ϵ * σ * (Ta + K0)^4
end


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
  φ₁ = 0.5 - 0.633 * χₗ - 0.330 * χₗ^2
  φ₂ = 0.877 * (1 - 2 * φ₁)
  # This equation is restricted to 0:4 < χℓ < 0:6
  φ₁ + φ₂ * cosZ # Bonan 2019, Eq. 14.31
end


function cal_τb_cum(dLAI::Vector{T}, Kb::T, Ω::T) where {T<:Real}
  n = length(dLAI)
  τb_cum = ones(T, n) # 第一层是地表
  # τb_cum[end] = 1.0 # 确认无误, 第N层必须是1
  ∑LAI = T(0.0)
  for i in n:-1:2
    ∑LAI += dLAI[i]
    τb_cum[i-1] = exp(-Kb * ∑LAI * Ω) # 确认无误
  end
  return τb_cum
end


# τ for Beam Radiation
function cal_τb(dLAI::Vector{T}, Kb::T, Ω::T) where {T<:Real}
  exp.(-Kb .* dLAI .* Ω)
end


# τ for Diffuse Radiation
function cal_τd(dLAI::Vector{T}, χₗ::T, Ω::T) where {T<:Real}
  φ₁ = 0.5 - 0.633 * χₗ - 0.330 * χₗ^2
  φ₂ = 0.877 * (1 - 2 * φ₁)

  τd = zeros(length(dLAI))
  for j in 1:9
    Zᵢ = (5 + (j - 1) * 10) |> deg2rad # 天顶角
    G_μ = φ₁ + φ₂ * cos(Zᵢ)
    τd .+= exp.(-G_μ / cos(Zᵢ) .* dLAI .* Ω) .* sin(Zᵢ) .* cos(Zᵢ)
  end
  τd .*= 2deg2rad(10)
  return τd
end


"""
阳叶比例, Kb = L_h / L, 将L转换为有效拦截面积

- `dLAI`: 0-N; 0对应地表
"""
function cal_fsun(dLAI::Vector{T}, Kb::T, Ω::T) where {T<:Real}
  _dLAI = @view dLAI[2:end]
  LAI = sum(_dLAI)
  sumlai = vcat(NaN, LAI .- cumsum(_dLAI) .+ _dLAI ./ 2)
  f_sun = @. Ω * exp(-Kb * Ω * sumlai) # 阳叶的比例, 光线被拦截的比例
  return f_sun
end
# laisun = (1 - exp(-Kb * LAI * Ω)) / Kb # 这可能是做了一个积分
# laisha = lai - laisun


function build_tridiag_shortwave!(
  tridiag::TriDiagonal{T}, nlayer::Int,
  Rs_dir::T, Rs_dif::T,
  τd::Vector{T}, τb::Vector{T}, τb_cum::Vector{T}, τl::T,
  ρ::T, ρ_soil_dir::T, ρ_soil_dif::T
) where {T}
  (; a, b, c, d) = tridiag # n = nlayer * 2 + 2  

  # Soil: upward
  m = 1
  iv = 1
  a[m] = T(0)
  b[m] = T(1)
  c[m] = -ρ_soil_dif
  d[m] = Rs_dir * τb_cum[m] * ρ_soil_dir

  # Soil: downward
  refld = (1 - τd[iv+1]) * ρ
  trand = (1 - τd[iv+1]) * τl + τd[iv+1]
  aiv = refld - trand^2 / refld
  biv = trand / refld

  m = 2
  a[m] = -aiv
  b[m] = T(1)
  c[m] = -biv
  d[m] = Rs_dir * τb_cum[iv+1] * (1 - τb[iv+1]) * (τl - ρ * biv)

  # Leaf layers
  @inbounds for iv in 2:nlayer
    # Upward
    refld = (1 - τd[iv]) * ρ
    trand = (1 - τd[iv]) * τl + τd[iv]
    fiv = refld - trand^2 / refld
    eiv = trand / refld

    m += 1
    a[m] = -eiv
    b[m] = T(1)
    c[m] = -fiv
    d[m] = Rs_dir * τb_cum[iv] * (1 - τb[iv]) * (ρ - τl * eiv)

    # Downward
    refld = (1 - τd[iv+1]) * ρ
    trand = (1 - τd[iv+1]) * τl + τd[iv+1]
    aiv = refld - trand^2 / refld
    biv = trand / refld

    m += 1
    a[m] = -aiv
    b[m] = T(1)
    c[m] = -biv
    d[m] = Rs_dir * τb_cum[iv+1] * (1 - τb[iv+1]) * (τl - ρ * biv)
  end

  # Top canopy layer: upward
  iv = nlayer + 1
  refld = (1 - τd[iv]) * ρ
  trand = (1 - τd[iv]) * τl + τd[iv]
  fiv = refld - trand^2 / refld
  eiv = trand / refld

  m += 1
  a[m] = -eiv
  b[m] = T(1)
  c[m] = -fiv
  d[m] = Rs_dir * τb_cum[iv] * (1 - τb[iv]) * (ρ - τl * eiv)

  # Top canopy layer: downward
  m += 1
  a[m] = T(0)
  b[m] = T(1)
  c[m] = T(0)
  d[m] = Rs_dif
end

function unpack_U_shortwave!(u::Vector{T}, I_up::Vector{T}, I_dn::Vector{T};
  nlayer::Int) where {T}
  I_up[1] = u[1] # I↑₀, I↓₀, ..., I↑ₙ, I↓ₙ
  I_dn[1] = u[2]

  # Leaf layer fluxes
  @inbounds for iv in 2:(nlayer+1)
    i1 = (iv - 1) * 2 + 1
    i2 = (iv - 1) * 2 + 2
    I_up[iv] = u[i1]
    I_dn[iv] = u[i2]
  end
end
