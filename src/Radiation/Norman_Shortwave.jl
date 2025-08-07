export Norman_Shortwave


@with_kw struct CanopyLeaves{FT}
  nlayer::Int = 2
  LAI::FT = 2.0
  Ω::FT = 1.0
  dLAI::Vector{FT} = ones(FT, nlayer) .* (LAI / nlayer)

  # 从下标从地面开始，从下至上
  LAI_sunlit::Vector{FT} = zeros(FT, nlayer)
  LAI_shaded::Vector{FT} = zeros(FT, nlayer)

  τ_d::Vector{FT} = zeros(FT, nlayer) # diffuse
  τ_b::Vector{FT} = zeros(FT, nlayer) # beam 

  Rs_up::Vector{FT} = zeros(FT, nlayer) # [I↑₁, I↑₂, ..., I↑ₙ]
  Rs_dn::Vector{FT} = zeros(FT, nlayer) # [I↓₁, I↓₂, ..., I↓ₙ]

  tridiag::TriDiagonal{FT} = TriDiagonal{FT}(; n=(nlayer + 1) * 2)
end


"""
    Norman_Shortwave(
        dLAI::Vector{T},
        Rs_dir=1000., Rs_dif=200.,
        τl::T=0.05, ρ::T=0.1, ρ_soil_dir::T=0.1, ρ_soil_dif::T=0.1,
        cosZ::T=0.88, χₗ::T=0.1, Ω::T=0.8) where {T<:Real}

# Arguments
- dLAI: 1-N, 注意是从下之上

# 遗留问题
1. 未考虑近红外和可见光的区别
2. 已知有叶片倾角数据，如何设定χₗ
3. 土壤对dir和dif的反射率
4. 叶片透射率
5. 找一个有多层Rs观测的站点，测试辐射传输
"""
function Norman_Shortwave(
  dLAI::Vector{T},
  Rs_dir=1000., Rs_dif=200.,
  τl=0.05, ρ=0.1, ρ_soil_dir=0.1, ρ_soil_dif=0.1,
  cosZ=0.88, χₗ=0.1, Ω=0.8) where {T<:Real}

  ϵ = 1 - (ρ + τl)
  nlayer = length(dLAI)
  dLAI = vcat(NaN, dLAI) # 增加地表

  Kb = min(cal_G(χₗ, cosZ) / cosZ, 20.0) # L_H/L
  f_sun = cal_fsun(dLAI, Kb, Ω)
  f_sha = 1 .- f_sun

  τb_cum = cal_τb_cum(dLAI, Kb, Ω)
  τb = cal_τb(dLAI, Kb, Ω)
  τd = cal_τd(dLAI, χₗ, Ω)

  tridiag = TriDiagonal{T}(; n=2 * (nlayer + 1))
  build_tridiag_shortwave!(
    tridiag, nlayer, Rs_dir, Rs_dif,
    τd, τb, τb_cum, τl,
    ρ, ρ_soil_dir, ρ_soil_dif
  )
  tridiagonal_solver!(tridiag)

  I_up = zeros(nlayer + 1) # Upward diffuse solar flux above layer
  I_dn = zeros(nlayer + 1) # Downward diffuse solar flux onto layer
  unpack_U_shortwave!(tridiag.u, I_up, I_dn; nlayer)

  # Absorbed direct beam and diffuse for each leaf layer and sum
  # for all leaf layers
  Rs_veg = T(0)
  Rs_sun = zeros(nlayer + 1)
  Rs_sha = zeros(nlayer + 1)

  for iv in 2:(nlayer+1)
    # Per unit ground area (W/m2 ground)
    direct = Rs_dir * τb_cum[iv] * (1 - τb[iv]) * ϵ
    diffuse = (I_dn[iv] + I_up[iv-1]) * (1 - τd[iv]) * ϵ

    # Absorbed solar radiation for shaded and sunlit portions of leaf layer
    # per unit ground area (W/m2 ground)
    sun = diffuse * f_sun[iv] + direct
    shade = diffuse * f_sha[iv]

    # Convert to per unit sunlit and shaded leaf area (W/m2 leaf)
    Rs_sun[iv] = sun / (f_sun[iv] * dLAI[iv])
    Rs_sha[iv] = shade / (f_sha[iv] * dLAI[iv])

    # Sum fluxes over all leaf layers
    Rs_veg = Rs_veg + (direct + diffuse) # note: 这是单位地面，和单位lai不同
  end
  # Rs_veg_sun = sum(Rs_sun)
  # Rs_veg_sha = sum(Rs_sha)
  direct = Rs_dir * τb_cum[1] * (1 - ρ_soil_dir)
  diffuse = I_dn[1] * (1 - ρ_soil_dif)
  Rs_gound = direct + diffuse

  # --- Conservation check
  # Total radiation balance: absorbed = incoming - outgoing
  incoming = Rs_dir + Rs_dif
  reflected = I_up[nlayer+1]
  sumabs = incoming - reflected

  err = sumabs - (Rs_veg + Rs_gound)
  if abs(err) > 1e-03
    println("err = %15.5f\n", err)
    error("NormanRadiation: Total solar conservation error")
  end

  (;
    PAR_sun=reverse(Rs_sun[2:(nlayer+1)]),
    PAR_sha=reverse(Rs_sha[2:(nlayer+1)]),
    frac_sha=reverse(f_sha[2:(nlayer+1)]),
    frac_sun=reverse(f_sun[2:(nlayer+1)])
  )
end
