import FieldMetadata: @bounds, bounds
import Ipaper: match2

abstract type AbstractETParam{FT<:AbstractFloat} end


"""
    struct Param_PMLV2{FT<:AbstractFloat} <: AbstractETParam{FT}

# Fields
$(TYPEDFIELDS)
"""
@bounds @with_kw mutable struct Param_PMLV2{FT<:Real} <: AbstractETParam{FT}
  "initial slope of the light response curve to assimilation rate, (i.e., quantum efficiency; `μmol CO2 [μmol PAR]⁻¹`)`"
  α::FT = 0.06 | (0.01, 0.10)
  "initial slope of the CO2 response curve to assimilation rate, (i.e., carboxylation efficiency; `μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹`)"
  η::FT = 0.04 | (0.01, 0.07)

  "stomatal conductance coefficient" # 气孔导度斜率参数
  g1::FT = 10.00 | (2.00, 100.00)

  "carbon saturated rate of photosynthesis at 25 °C, `μmol m⁻² s⁻¹`"
  Am_25::FT = 50.00 | (5.00, 120.00)

  "parameter to constrain `gc`, kPa"
  VPDmin::FT = 0.9 | (0.65, 1.5)
  "parameter to constrain `gc`, kPa"
  VPDmax::FT = 4.0 | (3.50, 6.5)

  "水汽压参数"  # leuning 2008
  D0::FT = 0.7 | (0.50, 2.0)
  "extinction coefficients for visible radiation" # 植被光合参数
  kQ::FT = 0.45 | (0.10, 1.0)
  "extinction coefficients for available energy"
  kA::FT = 0.70 | (0.50, 0.9)

  "Specific leaf storage, van Dijk, A.I.J.M, 2001, Eq2"
  S_sls::FT = 0.1 | (0.01, 1.0)
  "Canopy cover fraction related parameter"
  fER0::FT = 0.1 | (0.01, 0.5)

  "canopy height, `[m]`"
  hc::FT = 1.0 | (0.01, 20.0)

  "photoperiod constraint"
  d_pc::FT = 2.0 | (0.0, 5.0)

  ## 做出三套参数
  # LAIref::FT = 4.0   | (1.0, 6.0)      # 
  # frame::Integer = 10.0  | (6.0, 14.0) # 8-day moving window
  _hc::Vector{FT} = [1.0, 1.0, 1.0] | ([0.1, 0.1, 0.1], [20.0, 20.0, 20.0]) #(0.01, 20.0)
  _η::Vector{FT} = [0.04, 0.04, 0.04] | ([0.01, 0.01, 0.01], [0.07, 0.07, 0.07])
  _α::Vector{FT} = [0.06, 0.06, 0.06] | ([0.01, 0.01, 0.01], [0.10, 0.10])
  _g1::Vector{FT} = [10.0, 10.0, 10.0] | ([2.0, 2.0, 2.0], [100.0, 100.0, 100.0])
  _Am_25::Vector{FT} = [50.0, 50.0, 50.0] | ([5.0, 5.0, 5.0], [120.0, 120.0, 120.0])
end

## 划分为两种还是三种
ParNames = fieldnames(Param_PMLV2) |> collect

function Base.collect(par::AbstractETParam)
  [getfield(par, f) for f in fieldnames(typeof(par))]
end

function get_bounds(parNames::Vector{Symbol}=ParNames)
  inds = match2(parNames, ParNames).I_y # 选择的参数
  bs = bounds(Param_PMLV2)[inds]
  lower = vcat(map(x -> x[1], bs)...)
  upper = vcat(map(x -> x[2], bs)...)
  lower, upper
end


# theta2par(theta) = Param_PMLV2(theta...)
function theta2par!(theta::Vector, par::AbstractETParam, parNames::Vector{Symbol})
  k = 0
  for key in parNames
    x = getfield(par, key)
    n = length(x)
    inds = n == 1 ? k+1 : k+1:k+n
    k += n
    setfield!(par, key, theta[inds])
  end
  return par
end

function theta2par(theta::Vector, parNames::Vector{Symbol})
  theta2par!(theta, deepcopy(par0), parNames)
end

function select_param(par::AbstractETParam, parNames)
  vcat([getfield(par, f) for f in parNames]...)
end

# canopy height, in the order of `IGBP006` code
hc_raw = [10, 10, 10, 10, 10, 1, 1, 5, 5, 0.2, 1, 0.5, 10, 1, 0.01, 0.05, 0.01]

FT = Float64
par0 = Param_PMLV2()
theta0 = collect(par0)

function theta2param(params::Vector{Vector{T}}, IGBPs) where {T<:Real}
  parNames = fieldnames(Param_PMLV2) |> collect
  mat_param = map(collect, params) |> x -> cat(x..., dims=2) |> transpose
  d_param = DataFrame(mat_param, parNames)
  d_param.IGBP = IGBPs
  d_param[:, Cols(:IGBP, 1:end)]
end


export bounds, get_bounds, select_param
export Param_PMLV2, ParNames
export theta0, par0, hc_raw
export param_PML, theta2par, theta2param
