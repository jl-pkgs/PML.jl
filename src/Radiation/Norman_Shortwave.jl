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
- dLAI: 注意是从下之上
"""
function Norman_Shortwave(
  dLAI::Vector{T},
  Rs_dir=1000., Rs_dif=200.,
  ρ=0.1, τ_l=0.05,
  ρ_soil_dir=0.1, ρ_soil_dif=0.1,
  cosZ=0.88, χₗ=0.1, Ω=0.8) where {T<:Real}

  ϵ = 1 - (ρ + τ_l)
  nlayer = length(dLAI)

  lai = sum(dLAI)
  sumlai = vcat(NaN, lai .- cumsum(dLAI) .+ dLAI ./ 2)
  dLAI = vcat(NaN, dLAI) # 增加地表

  Kb = min(cal_G(χₗ, cosZ) / cosZ, 20.0)
  f_sun = Ω * exp.(-Kb * sumlai * Ω)
  f_sha = 1 .- f_sun
  # laisun = (1 - exp(-Kb * lai * Ω)) / Kb
  # laisha = lai - laisun

  τb_cum = cal_τb_cum(Kb, Ω, dLAI)
  τb = cal_τb(Kb, Ω, dLAI)
  τd = cal_τd(χₗ, Ω, dLAI)

  tridiag = TriDiagonal{T}(; n=2 * (nlayer + 1))

  build_tridiag_shortwave!(
    tridiag, nlayer, τd, τb, τb_cum, 
    Rs_dir, Rs_dif, 
    ρ, τ_l,
    ρ_soil_dir, ρ_soil_dif
  )
  tridiagonal_solver!(tridiag)

  I_up = zeros(nlayer + 1) # Upward diffuse solar flux above layer
  I_dn = zeros(nlayer + 1) # Downward diffuse solar flux onto layer
  unpack_U_shortwave!(tridiag.u, I_up, I_dn; nlayer)

  # --- Compute flux densities
  # Absorbed direct beam and diffuse for ground (soil)
  iv = 1
  direct = Rs_dir * τb_cum[iv] * (1 - ρ_soil_dir)
  diffuse = I_dn[iv] * (1 - ρ_soil_dif)
  Rs_gound = direct + diffuse

  # Absorbed direct beam and diffuse for each leaf layer and sum
  # for all leaf layers
  Rs_veg = T(0)
  Rs_veg_sun = T(0)
  Rs_veg_sha = T(0)
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
    Rs_veg = Rs_veg + (direct + diffuse)
    Rs_veg_sun = Rs_veg_sun + sun
    Rs_veg_sha = Rs_veg_sha + shade
  end

  # --- Albedo
  incoming = Rs_dir + Rs_dif
  reflected = I_up[nlayer+1]
  albcan = incoming > 0 ? reflected / incoming : 0

  # --- Conservation check
  # Total radiation balance: absorbed = incoming - outgoing
  suminc = Rs_dir + Rs_dif
  sumref = albcan * (Rs_dir + Rs_dif)
  sumabs = suminc - sumref

  err = sumabs - (Rs_veg + Rs_gound)
  if abs(err) > 1e-03
    println("err = %15.5f\n", err)
    error("NormanRadiation: Total solar conservation error")
  end

  # Sunlit and shaded absorption
  err = (Rs_veg_sun + Rs_veg_sha) - Rs_veg
  if abs(err) > 1e-03
    println("err = %15.5f\n", err)
    error("NormanRadiation: Sunlit/shade solar conservation error")
  end

  (;
    PAR_sun=reverse(Rs_sun[2:(nlayer+1)]),
    PAR_sha=reverse(Rs_sha[2:(nlayer+1)]),
    frac_sha=reverse(f_sha[2:(nlayer+1)]),
    frac_sun=reverse(f_sun[2:(nlayer+1)])
  )
end
