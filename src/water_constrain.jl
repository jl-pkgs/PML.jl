@with_kw struct β_GPP_Zhang2019{FT} <: AbstractWaterConsGPPModel{FT}
  ## water constraint
  "parameter to constrain `gc`, kPa"
  VPDmin::Param = Param(0.9, bounds=(0.65, 1.5))

  "parameter to constrain `gc`, kPa"
  VPDmax::Param = Param(4.0, bounds=(3.50, 6.5))
end


"Piecewise function by Yongqiang, Dongdong and GanRong, 2019"
function β_GPP(VPD::T, par::β_GPP_Zhang2019{T}) where {T<:Real}
  (; VPDmin, VPDmax) = par
  if (VPD > VPDmax)
    T(0.0)
  elseif VPD < VPDmin
    T(1.0)
  else
    (VPDmax - VPD) / (VPDmax - VPDmin)
  end
end


function f_VPD_Zhang2019(VPD::T, VPDmin::T, VPDmax::T)::T where {T<:Real}
  if (VPD > VPDmax)
    T(0.0)
  elseif VPD < VPDmin
    T(1.0)
  else
    (VPDmax - VPD) / (VPDmax - VPDmin)
  end
end
