"""
    PMLV2 (Penman–Monteith–Leuning Version 2) Evapotranspiration model

# Arguments

- `Prcp` : mm/d
- `Tavg` : degC
- `Rs`   : W m-2
- `Rn`   : W m-2
- `VPD`  : W m-2
- `U2`   : m/s
- `LAI`  : m2 m-2
- `Pa`   : kPa
- `Ca`   : ppm, default `380`
- `Ω`    : clamping index, default is `1.0`

# Examples
```julia
```

# References
1. Gan Rong, 2018, Ecohydrology
2. Zhang Yongqiang, 2019, RSE
3. Kong Dongdong, 2019, ISPRS
"""
function PMLV2(Prcp::T, Tavg::T, Rs::T, Rn::T, VPD::T, U2::T, LAI::T,
  Pa=atm, 
  Ca=380.0, PC=1.0,
  Ω::T=T(1.0);
  # leaf::AbstractLeaf, 
  par::Param_PMLV2=Param_PMLV2(),
  r::Union{Nothing,interm_PML}=nothing) where {T<:Real}
  r === nothing && (r = interm_PML{T}())

  # D0 = 0.7
  # kQ = 0.6 # extinction coefficients for visible radiation
  # kA = 0.7 # extinction coefficient for available energy
  λ, Δ, γ, r.Eeq = ET0_eq(Rn, Tavg, Pa)
  ϵ = Δ / γ

  ### CARBON MODULE: PHOTOSYNTHESIS --------------------------------------------
  r.GPP, r.Gc_w = photosynthesis(Tavg, Rs, VPD, LAI, Pa, Ca, PC; par)

  ### WATER MODULE: ------------------------------------------------------------
  ## Intercepted Evaporation (Ei)
  r.Ei = cal_Ei_Dijk2021(Prcp, LAI, par)
  r.Pi = Prcp - r.Ei # 应扣除这一部分消耗的能量

  r.Ga = aerodynamic_conductance(U2, par.hc) # Leuning, 2008, Eq.13, doi:10.1029/2007WR006562

  # Transpiration from plant cause by radiation water transfer
  Tou = exp(-par.kA * LAI)
  LEcr = ϵ * Rn * (1 - Tou) / (ϵ + 1 + r.Ga / r.Gc_w) # W m-2

  # Transpiration from plant cause by aerodynamic water transfer
  ρa = cal_rho_a(Tavg, Pa)
  # ρa = cal_rho_a(Tavg, q, Pa)
  LEca = (ρa * Cp * 1e6 * r.Ga * VPD / γ) / (ϵ + 1 + r.Ga / r.Gc_w) # W m-2, `Cp*1e6`: [J kg-1 °C-1]
  # [kg m-3] [J kg-1 K-1] [m s-1] [kPa] / [kPa K-1] = [W m-2]

  r.Ecr = W2mm(LEcr; λ) # [W m-2] change to [mm d-1]
  r.Eca = W2mm(LEca; λ) # [W m-2] change to [mm d-1]
  r.Ec = r.Ecr + r.Eca

  ## TODO: 补充冰面蒸发的计算
  Evp::T = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ
  r.ET_water = r.Eeq + Evp

  r.Es_eq = r.Eeq * Tou # Soil evaporation at equilibrium, mm d-1
  # r.Es_eq = r.Eeq * exp(-0.9 * LAI) # Ka = 0.9, insensitive param
  r
  # GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, Eeq, ET_water, Ga, Gc_w
end



"""
    PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca; par=param0, frame=3)

# Notes
一个站点的计算。注意，不同植被类型，参数不同。

# Arguments
- `frame`: in 8-days
"""
function PMLV2(Prcp::V, Tavg::V, Rs::V, Rn::V,
  VPD::V, U2::V, LAI::V,
  Pa::V,
  Ca::V,
  PC::Union{T,V} = T(1.0);
  par::Param_PMLV2=Param_PMLV2(), frame=3,
  res::Union{Nothing,output_PML}=nothing) where {T<:Real,V<:AbstractVector{T}}

  n = length(Prcp)
  fields = fieldnames(interm_PML)[2:end-2] # skip ET, fval_soil and Es
  r = interm_PML{T}()
  res === nothing && (res = output_PML{T}(; n))

  # isvec_Ca = isa(Ca, AbstractVector)
  isvec_PC = isa(PC, AbstractVector)
  @inbounds for t = 1:n
    # _Ca = isvec_Ca ? Ca[t] : Ca
    _PC = isvec_PC ? PC[t] : PC

    PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t],
      Ca[t], _PC; par, r)
    res[t, fields] = r
  end

  res.fval_soil = movmean2(Prcp, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.fval_soil, 0, 1)

  res.Es .= res.fval_soil .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end


"""
# Arguments
- `kw`: named keyword arguments
  + `r`: `interm_PML`
"""
function PMLV2(d::AbstractDataFrame; par::Param_PMLV2=Param_PMLV2(), kw...)
  (; Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca) = d
  PC = "PC" ∈ names(d) ? d.PC : 1.0
  
  PMLV2(Prcp, Tavg, Rs, Rn,
    VPD, U2, LAI,
    Pa, Ca, PC; par, kw...) |> to_df
end

# 相同植被类型多个站点一起的计算
function PMLV2_sites(df::AbstractDataFrame; par::Param_PMLV2=Param_PMLV2(), kw...)
  "site" ∉ names(df) && return PMLV2(df; par, kw...)

  sites = df.site
  grps = unique(sites)

  # ! 这里有优化空间，但程序会写的很复杂
  res = []
  for grp in grps
    inds = sites .== grp
    d = df[inds, :]
    r = PMLV2(d; par, kw...)
    push!(res, r)
  end
  vcat(res...)
end
