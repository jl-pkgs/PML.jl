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

