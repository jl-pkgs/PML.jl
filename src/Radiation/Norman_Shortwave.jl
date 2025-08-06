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


function Norman_Shortwave(
  dLAI,
  PARdir=1000, PARdif=200,
  ρ=0.1, τ_l=0.05,
  ρ_soil_dir=0.1, ρ_soil_dif=0.1,
  cosz=0.88, chil=0.1, Ω=0.8)

  nlayer = length(dLAI)
  dLAI = reverse(dLAI) # TODO: 为何要反转?

  lai = sum(dLAI)
  sumlai = vcat(NaN, lai .- cumsum(dLAI) .+ dLAI ./ 2)
  dLAI = vcat(NaN, dLAI)

  if chil > 0.6 || chil < -0.4
    println("Chil is not inside the interval -0.4, 0.6 and was changed")
  end

  chil = clamp(chil, -0.4, 0.6)
  phi1 = 0.5 - 0.633 * chil - 0.330 * chil^2
  phi2 = 0.877 * (1 - 2 * phi1)

  gdir = phi1 + phi2 * cosz

  Kb = gdir / cosz
  Kb = min(Kb, 20)

  fracsun = Ω * exp.(-Kb * sumlai * Ω)
  fracsha = 1 .- fracsun
  # laisun = (1 - exp(-Kb * lai * Ω)) / Kb
  # laisha = lai - laisun

  τ_b = exp.(-Kb * dLAI * Ω)
  τ_d = zeros(length(dLAI))
  for j in 1:9
    angle = (5 + (j - 1) * 10) * pi / 180
    G_μ = phi1 + phi2 * cos(angle)
    τ_d = τ_d .+ exp.(-G_μ / cos(angle) * dLAI * Ω) .* sin(angle) .* cos(angle)
  end
  τ_d = τ_d .* 2 .* (10 * pi / 180)

  τ_bcum = vcat(fill(NaN, nlayer + 1))
  cumlai = 0
  iv = nlayer + 1
  τ_bcum[iv] = 1
  for iv in (nlayer+1):-1:2
    cumlai += dLAI[iv]
    τ_bcum[iv-1] = exp(-Kb * cumlai * Ω)
  end
  # println("Radiation model for a total LAI of ", lai)

  swup = zeros(nlayer + 1)
  swdn = zeros(nlayer + 1)
  a = zeros(nlayer * 2 + 2)
  b = zeros(nlayer * 2 + 2)
  c = zeros(nlayer * 2 + 2)
  d = zeros(nlayer * 2 + 2)

  ϵ = 1 - (ρ + τ_l)
  m = 1
  iv = 1
  a[m] = 0
  b[m] = 1
  c[m] = -ρ_soil_dif
  d[m] = PARdir * τ_bcum[m] * ρ_soil_dir

  # Soil: downward flux
  refld = (1 - τ_d[iv+1]) * ρ
  trand = (1 - τ_d[iv+1]) * τ_l + τ_d[iv+1]
  aiv = refld - trand * trand / refld
  biv = trand / refld

  m = 2
  a[m] = -aiv
  b[m] = 1
  c[m] = -biv
  d[m] = PARdir * τ_bcum[iv+1] * (1 - τ_b[iv+1]) * (τ_l - ρ * biv)

  # 这里有优化的空间
  # Leaf layers, excluding top layer
  for iv in 2:nlayer
    # Upward flux
    refld = (1 - τ_d[iv]) * ρ
    trand = (1 - τ_d[iv]) * τ_l + τ_d[iv]
    fiv = refld - trand * trand / refld
    eiv = trand / refld

    m += 1
    a[m] = -eiv
    b[m] = 1
    c[m] = -fiv
    d[m] = PARdir * τ_bcum[iv] * (1 - τ_b[iv]) * (ρ - τ_l * eiv)

    # Downward flux
    refld = (1 - τ_d[iv+1]) * ρ
    trand = (1 - τ_d[iv+1]) * τ_l + τ_d[iv+1]
    aiv = refld - trand * trand / refld
    biv = trand / refld

    m += 1
    a[m] = -aiv
    b[m] = 1
    c[m] = -biv
    d[m] = PARdir * τ_bcum[iv+1] * (1 - τ_b[iv+1]) * (τ_l - ρ * biv)
  end

  # Top canopy layer: upward flux
  iv = nlayer + 1
  refld = (1 - τ_d[iv]) * ρ
  trand = (1 - τ_d[iv]) * τ_l + τ_d[iv]
  fiv = refld - trand * trand / refld
  eiv = trand / refld

  m += 1
  a[m] = -eiv
  b[m] = 1
  c[m] = -fiv
  d[m] = PARdir * τ_bcum[iv] * (1 - τ_b[iv]) * (ρ - τ_l * eiv)

  # Top canopy layer: downward flux
  m += 1
  a[m] = 0
  b[m] = 1
  c[m] = 0
  d[m] = PARdif
  u = tridiagonal_solver(a, b, c, d) # Solve tridiagonal equations for fluxes

  # Now copy the solution (u) to the upward (swup) and downward (swdn) fluxes for each layer
  # swup - Upward diffuse solar flux above layer
  # swdn - Downward diffuse solar flux onto layer
  # Soil fluxes
  iv = 1
  m = 1
  swup[iv] = u[m]
  m += 1
  swdn[iv] = u[m]

  # Leaf layer fluxes
  for iv in 2:(nlayer+1)
    i1 = (iv - 1) * 2 + 1
    i2 = (iv - 1) * 2 + 2
    swup[iv] = u[i1]
    swdn[iv] = u[i2]
  end

  # --- Compute flux densities
  # Absorbed direct beam and diffuse for ground (soil)
  iv = 1
  direct = PARdir * τ_bcum[iv] * (1 - ρ_soil_dir)
  diffuse = swdn[iv] * (1 - ρ_soil_dif)
  swsoi = direct + diffuse

  # Absorbed direct beam and diffuse for each leaf layer and sum
  # for all leaf layers
  swveg = 0
  swvegsun = 0
  swvegsha = 0
  Rs_sun = zeros(nlayer + 1)
  Rs_sha = zeros(nlayer + 1)

  for iv in 2:(nlayer+1)
    # Per unit ground area (W/m2 ground)
    direct = PARdir * τ_bcum[iv] * (1 - τ_b[iv]) * ϵ
    diffuse = (swdn[iv] + swup[iv-1]) * (1 - τ_d[iv]) * ϵ

    # Absorbed solar radiation for shaded and sunlit portions of leaf layer
    # per unit ground area (W/m2 ground)
    sun = diffuse * fracsun[iv] + direct
    shade = diffuse * fracsha[iv]

    # Convert to per unit sunlit and shaded leaf area (W/m2 leaf)
    Rs_sun[iv] = sun / (fracsun[iv] * dLAI[iv])
    Rs_sha[iv] = shade / (fracsha[iv] * dLAI[iv])

    # Sum fluxes over all leaf layers
    swveg = swveg + (direct + diffuse)
    swvegsun = swvegsun + sun
    swvegsha = swvegsha + shade
  end

  # --- Albedo
  incoming = PARdir + PARdif
  reflected = swup[nlayer+1]
  albcan = incoming > 0 ? reflected / incoming : 0

  # --- Conservation check
  # Total radiation balance: absorbed = incoming - outgoing
  suminc = PARdir + PARdif
  sumref = albcan * (PARdir + PARdif)
  sumabs = suminc - sumref

  err = sumabs - (swveg + swsoi)
  if abs(err) > 1e-03
    println("err = %15.5f\n", err)
    error("NormanRadiation: Total solar conservation error")
  end

  # Sunlit and shaded absorption
  err = (swvegsun + swvegsha) - swveg
  if abs(err) > 1e-03
    println("err = %15.5f\n", err)
    error("NormanRadiation: Sunlit/shade solar conservation error")
  end

  (;
    PAR_sun=reverse(Rs_sun[2:(nlayer+1)]),
    PAR_sha=reverse(Rs_sha[2:(nlayer+1)]),
    frac_sha=reverse(fracsha[2:(nlayer+1)]),
    frac_sun=reverse(fracsun[2:(nlayer+1)])
  )
end



# function update_ef2(atri, btri, ctri, dtri, ir_source, td, iv, m)
#   refld = (1 - td[iv]) * ρ
#   trand = (1 - td[iv]) * τ + td[iv]
#   fiv = refld - trand * trand / refld
#   eiv = trand / refld

#   atri[m] = -eiv
#   btri[m] = 1
#   ctri[m] = -fiv
#   dtri[m] = (1 - eiv) * ir_source[iv]
#   nothing
#   # fiv, eiv
# end
