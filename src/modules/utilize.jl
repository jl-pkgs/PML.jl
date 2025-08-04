# 最大羧化能力温度调节函数
# V_m = Vm_25 * T_adjust_Vm25
# `T` in -100:100, `T_adjust_Vm25` always < 0.92
function T_adjust_Vm25(Tavg::T)::T where {T<:Real}
  a = 0.031
  b = 0.115
  exp(a * (Tavg - 25.0)) / (1.0 + exp(b * (Tavg - 41.0))) # Gan2018, Eq. A5
end


const R_gas = 8.3144621 # J/(mol K)

"""
    mol2m(x, Tavg, Pa=atm)
    mol2m_rong2018(x, Tavg, Pa=atm)

Convert from umol m-2 s-1 to m s-1, g = g_m * mol2m(Tavg)

## Arguments
- `Tavg`: degC
- `Pa`: kPa

# Reference
1. Monteith, 2013, Principles of Environmental Physics, Eq. 3.14
1. CLM5 TechNote, Eq. 9.1
"""
mol2m(x::T, Tavg::T, Pa::T=T(atm)) where {T} = x * T(R_gas) * (Tavg + K0) / (Pa * 1e3) # the last in [umol to ]

mol2m_rong2018(x::T, Tavg::T, Pa::T=T(atm)) where {T} =
  x * 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3))

# [umol m-2 s-1] to [g C m-2 d-1]
umol2gC(mol::T) where {T} = mol * 86400 / 10^6 * 12

gC2umol(x::T) where {T} = x / (12 * 0.0864)

export mol2m, mol2m_rong2018, umol2gC, gC2umol
