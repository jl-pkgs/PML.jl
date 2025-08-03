abstract type AbstractβGPPModel end


@with_kw struct β_GPP_Zhang2019{FT} <: AbstractβGPPModel{FT}
  ## water constraint
  "parameter to constrain `gc`, kPa"
  VPDmin::FT = Param(0.9, bounds=(0.65, 1.5))

  "parameter to constrain `gc`, kPa"
  VPDmax::FT = Param(4.0, bounds=(3.50, 6.5))
end


function β_GPP(VPD::T, par::β_GPP_Zhang2019)
  (; VPDmin, VPDmax) = par
  if (VPD > VPDmax)
    T(0.0)
  elseif VPD < VPDmin
    T(1.0)
  else
    (VPDmax - VPD) / (VPDmax - VPDmin)
  end
end


# Piecewise function by Yongqiang and GanRong, 2019
function f_VPD_Zhang2019(VPD::T, VPDmin::T, VPDmax::T)::T where {T<:Real}
  if (VPD > VPDmax)
    T(0.0)
  elseif VPD < VPDmin
    T(1.0)
  else
    (VPDmax - VPD) / (VPDmax - VPDmin)
  end
end
