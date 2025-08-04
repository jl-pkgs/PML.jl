export β_GPP_Zhang2019


@bounds @units @with_kw mutable struct β_GPP_Zhang2019{FT} <: AbstractWaterConsGPPModel{FT}
  ## water constraint
  "parameter to constrain `gc`"
  VPDmin::FT = 0.9 | (0.65, 1.5) | "kPa"

  "parameter to constrain `gc`"
  VPDmax::FT = 4.0 | (3.50, 6.5) | "kPa"
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
