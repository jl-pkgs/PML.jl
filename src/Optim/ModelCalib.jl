export ModelCalib, ModelCalib_IGBPs
export fread, fwrite, melt_list

using RTableTools
import Ipaper: par_map
include("DataType.jl")
include("model_gof.jl")
# include("PMLV2_sites.jl")


function ModelCalib_IGBPs(data::AbstractDataFrame; of_gof=:KGE, maxn=2500)
  vars = ["IGBPname", "IGBPcode", "site", "date", "GPP_obs", "ET_obs",
    "Prcp", "Tavg", "U2", "Rn", "Rs", "VPD", "LAI", "Pa", "Ca"]
  IGBPs = unique(data.IGBPname) |> sort

  @time params = par_map(IGBP -> begin
      df = data[data.IGBP.==IGBP, vars]
      IGBPcode = df.IGBPcode[1]
      theta, _, _ = ModelCalib(df, par0; IGBPcode, of_gof, maxn, verbose=false)
      theta
    end, IGBPs; progress=false)

  # printstyled("[i=$i, IGBP = $IGBP] \n", bold=true, color=:green, underline=false)
  res = map(i -> begin
      IGBP = IGBPs[i]
      theta = params[i]
      df = data[data.IGBP.==IGBP, vars]
      model_goal(df, theta; verbose=true)

      par = theta2par(theta)
      r = PMLV2_sites(df; par)
      cbind(df[:, [:site, :date, :ET_obs, :GPP_obs]], r)
    end, eachindex(params))
  df_out = melt_list(res; IGBP=IGBPs)

  gof = model_gof(df_out)
  (; param=theta2param(params, IGBPs), gof, output=df_out)
end


## 一个站点的率定
"""
    ModelCalib(df::AbstractDataFrame, par0::AbstractETParam; 
        IGBPcode=nothing, maxn=2500, of_gof=:NSE, kw...)
"""
function ModelCalib(df::AbstractDataFrame, par0::AbstractETParam; IGBPcode=nothing, maxn=2500, of_gof=:NSE, kw...)
  parRanges = get_bounds(par0)
  lower = parRanges[:, 1]
  upper = parRanges[:, 2]
  theta0 = collect(par0) # must be vector

  sceua(theta -> -model_goal(df, theta; IGBPcode, of_gof, kw...),
    theta0, lower, upper; maxn, kw...) # theta, goal, flag
end
