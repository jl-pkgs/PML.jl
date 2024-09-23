# Piecewise function by Yongqiang and GanRong, 2019
function f_VPD_Zhang2019(VPD::T, par::Param_PMLV2)::T where {T<:Real}
  if (VPD > par.VPDmax)
    T(0.0)
  elseif VPD < par.VPDmin
    T(1.0)
  else
    (par.VPDmax - VPD) / (par.VPDmax - par.VPDmin)
  end
end
