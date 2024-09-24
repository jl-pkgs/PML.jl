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
