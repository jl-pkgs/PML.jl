import FieldMetadata: @bounds, bounds
import Ipaper: match2



## 做出三套参数
# LAIref::FT = 4.0   | (1.0, 6.0)      # 
# frame::Integer = 10.0  | (6.0, 14.0) # 8-day moving window

## 做出三套参数
# _hc::Vector{FT} = [1.0, 1.0, 1.0] | ([0.1, 0.1, 0.1], [20.0, 20.0, 20.0]) #(0.01, 20.0)
# _η::Vector{FT} = [0.04, 0.04, 0.04] | ([0.01, 0.01, 0.01], [0.07, 0.07, 0.07])
# _α::Vector{FT} = [0.06, 0.06, 0.06] | ([0.01, 0.01, 0.01], [0.10, 0.10, 0.10])
# _g1::Vector{FT} = [10.0, 10.0, 10.0] | ([2.0, 2.0, 2.0], [100.0, 100.0, 100.0])
# _VCmax25::Vector{FT} = [50.0, 50.0, 50.0] | ([5.0, 5.0, 5.0], [120.0, 120.0, 120.0])


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
export theta2par, theta2param
