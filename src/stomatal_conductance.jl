@bounds @with_kw mutable struct SM_Yu2004{FT} <: AbstractStomatalModel{FT}
  "stomatal conductance coefficient, `μmol m⁻² s⁻¹`" # 气孔导度斜率参数
  g1::FT = 10.00 | (2.00, 100.00)

  "水汽压参数"  # leuning 2008, Yu 2004
  D0::FT = 0.7 | (0.50, 2.0) # kPa
end


@bounds @with_kw struct SM_MedlynSM2011{FT} <: AbstractStomatalModel{FT}
  "stomatal conductance coefficient, `μmol m⁻² s⁻¹`" # 气孔导度斜率参数
  g0::FT = 100.00 | (50.0, 150.0)

  "stomatal conductance coefficient, `sqrt(kPa)`" # 气孔导度斜率参数
  g1::FT = 2.00 | (1.0, 6.0)
end



"""
- `An`: net assimilation rate, [umol m-2 s-1], `An = Ag - Rd`
- `An`: net assimilation rate, [umol m-2 s-1]
- `D` : air vapor pressure, [kPa]
- `Ca`: ~380 ppm, 380 [10^-6 mol mol-1]

> 1.6是`Gc` to `Gc_w`
"""
function stomatal_conductance(stomatal::SM_MedlynSM2011, Ag::T, Rd::T, D::T, Ca::T) where {T<:Real}
  (; g0, g1) = stomatal
  An = Ag - Rd
  gs = g0 + 1.6g1 * (1 + g1 / sqrt(D)) * An / Ca
  return gs
end


function stomatal_conductance(stomatal::SM_Yu2004, Ag::T, Rd::T, D::T, Ca::T) where {T<:Real}
  (; g1, D0) = stomatal
  gs = 1.6g1 * Ag / (Ca * (1 + D / D0)) # Yu 2004; Rong 2018 Eq. A3
  return gs
end
