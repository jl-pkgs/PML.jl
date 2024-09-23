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
