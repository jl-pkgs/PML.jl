abstract type AbstractStomataModel end


@with_kw struct SM_Yu2004{FT} <: AbstractStomataModel{FT}
  "stomatal conductance coefficient, `μmol m⁻² s⁻¹`" # 气孔导度斜率参数
  g1::FT = Param(10.00, bounds=(2.00, 100.00))

  "水汽压参数"  # leuning 2008, Yu 2004
  D0::FT = Param(0.7, bounds=(0.50, 2.0)) # kPa
end


@with_kw struct SM_MedlynSM2011{FT} <: AbstractStomataModel{FT}
  "stomatal conductance coefficient, `μmol m⁻² s⁻¹`" # 气孔导度斜率参数
  g0::FT = Param(100.00, bounds=(50.0, 150.0))

  "stomatal conductance coefficient, `sqrt(kPa)`" # 气孔导度斜率参数
  g1::FT = Param(2.00, bounds=(1.0, 6.0))
end



"""
- `An`: net assimilation rate, [umol m-2 s-1], `An = Ag - Rd`
- `An`: net assimilation rate, [umol m-2 s-1]
- `D` : air vapor pressure, [kPa]
- `Ca`: ~380 ppm, 380 [10^-6 mol mol-1]
"""
function stomatal_conductance(Ag::T, Rd::T, D::T, Ca::T, par::SM_MedlynSM2011) where {T<:Real}
  (; g0, g1) = par
  An = Ag - Rd
  g0 + 1.6g1 * (1 + g1 / sqrt(D)) * An / Ca
end


function stomatal_conductance(Ag::T, Rd::T, D::T, Ca::T, par::SM_Yu2004) where {T<:Real}
  (; g1, D0) = par
  g1 * Ag / (Ca * (1 + D / D0)) # Yu 2004; Rong 2018 Eq. A3
end
