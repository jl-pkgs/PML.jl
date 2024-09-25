import HydroTools: cal_slope

"""
    cal_Ei_Dijk2021(Prcp::T, LAI::T, par::Param_PMLV2) where {T<:Real}

# References
1. van Dijk, A.I.J.M, 2001, Eq2.
"""
function cal_Ei_Dijk2021(Prcp::T, LAI::T, par::Param_PMLV2) where {T<:Real}
  # two params in Ei
  # @unpack S_sls, fER0 = par
  LAIref = 5
  # van Dijk, A.I.J.M, 2001, Eq2.
  fveg = 1 - exp(-LAI / LAIref)  # Canopy cover fraction, Eq.1
  Sveg = par.S_sls * LAI # Specific leaf storage, Eq.2
  fER = par.fER0 * fveg # the value of 0.50 based on optimisation at Australian catchments
  Pwet = -log(1 - par.fER0) / par.fER0 * Sveg / fveg # -log(1 - fER /fveg),

  # Pwet[is.na(Pwet)] = 0; check, negative value, log will give error
  Ei = Prcp < Pwet ? fveg * Prcp : (fveg * Pwet + fER * (Prcp - Pwet))
  return Ei
end


"""
    aerodynamic_conductance(U2, hc)

# Arguments
- `U2`: wind speed at 2m
- `hc`: canopy height

# Return
- `Ga`: aerodynamic conductance in m/s
"""
function aerodynamic_conductance(U2::T, hc::T; Zob=15.0) where {T<:Real}
  kmar = 0.40        # von Karman's constant 0.40
  d = 0.64 * hc
  zom = 0.13 * hc
  zoh = 0.10 * zom
  uz = cal_Uz(U2, Zob) # convert from u2 to uz
  @fastmath Ga = uz * kmar^2 / (log((Zob - d) / zom) * log((Zob - d) / zoh)) # m s-1
  Ga
end


function ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...) where {T<:Real}
  # T = eltype(Rn)
  ϵ = 0.622 # ratio of molecular weight of water vapor to dry air
  Cp = 4.2 * 0.242 / 1000 # MJ kg-1 K-1
  
  λ::T = cal_lambda(Tair) # MJ kg-1
  Δ::T = cal_slope(Tair) # kPa degC-1
  γ::T = Cp * Pa / (ϵ * λ) # kPa degC-1
  Eeq::T = Δ / (Δ + γ) * Rn |> x -> W2mm(x; λ)
  λ, Δ, γ, Eeq
end

# lambda: [MJ kg-1]
W2mm(Ra; λ) = Ra * 86400 / 1e6 / λ

function cal_lambda(Tair::T) where {T<:Real}
  #  * u"MJ / kg"
  # 2.501 - 0.00237 * Tair # bolton 1980
  (2500.0 - Tair * 2.2) / 1000
end

cal_rho_a(Tair, Pa) = 3.486 * Pa / 1.01(Tair + 273.15) # kg/m3
# cal_rho_a(Tair, Pa) = 3.846 * Pa / (Tair + 273.15) # [kg/m3], error v2019

# # rho_a: kg m-3
# function cal_rho_a(Tair, q, Pa)
#   3.486 * Pa / cal_TvK(Tair, q) # FAO56, Eq. 3-5
# end

# cal_rho_a(Tair, Pa) = 3.486 * Pa / cal_TvK(Tair)
# # cal_rho_a(Tair, Pa) = 3.486 * Pa / cal_TvK(Tair)

# # FAO56, Eq. 3-7
# cal_TvK(Tair) = 1.01 * (Tair + 273)

# # https://github.com/CUG-hydro/class2022_CUG_HydroMet/blob/master/ch02_大气的基本特征.md
# # q ≈ ϵ*ea/Pa
# # q = ϵ*ea/(Pa - (1 - ϵ)*ea)
# function cal_TvK(Tair, q)
#   # ea / Pa = q/(ϵ + (1 - ϵ) * q)
#   (Tair + K0) * (1 + (1 - epsilon) * q / (ϵ + (1 - ϵ) * q))
# end

# # 这个是最精确的版本
# # FAO56, Eq. 3-6
# cal_TvK(Tair, ea, Pa) = (Tair + K0) * (1 + (1 - epsilon) * ea / Pa)

export cal_Ei_Dijk2021, aerodynamic_conductance
