export ModelCalib, ModelCalib_IGBPs
export model_goal
export model_gof

include("DataType.jl")

function ModelCalib_IGBPs(data::AbstractDataFrame;
  parNames=ParNames, of_gof=:KGE, maxn=2500)

  vars = ["IGBPname", "IGBPcode", "site", "date", "GPP_obs", "ET_obs",
    "Prcp", "Tavg", "U2", "Rn", "Rs", "VPD", "LAI", "Pa", "Ca"]
  IGBPs = unique(data.IGBPname) |> sort

  @time params = par_map(IGBP -> begin
      df = data[data.IGBP.==IGBP, vars]
      IGBPcode = df.IGBPcode[1]
      ModelCalib(df, par0, parNames; IGBPcode, of_gof, maxn, verbose=false)[1]
    end, IGBPs; progress=false)

  # printstyled("[i=$i, IGBP = $IGBP] \n", bold=true, color=:green, underline=false)
  res = map(i -> begin
      IGBP = IGBPs[i]
      theta = params[i]
      df = data[data.IGBP.==IGBP, vars]
      model_goal(df, theta, parNames; verbose=true) # print GOF

      par = theta2par(theta, parNames)
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
function ModelCalib(df::AbstractDataFrame, par0::AbstractETParam, parNames;
  IGBPcode=nothing, maxn=2500, of_gof=:NSE, kw...)
  lower, upper = get_bounds(parNames)
  x0 = select_param(par0, parNames)

  sceua(theta -> -model_goal(df, theta, parNames; IGBPcode, of_gof, kw...),
    x0, lower, upper; maxn, kw...) # theta, goal, flag
end


## gof of ET and GPP for one IGBP
# cal df_out by Forcing and Param
function model_goal(df, theta, parNames; IGBPcode=nothing, of_gof=:NSE, verbose=false)
  par = theta2par(theta, parNames)
  :hc ∉ parNames && !isnothing(IGBPcode) && (par.hc = hc_raw[IGBPcode])

  dobs = df[!, [:GPP_obs, :ET_obs]]
  dsim = PMLV2_sites(df; par)
  # dsim = PMLV2(df; par) # the exact MATLAB version

  ## 8-day, yearly, yearly_anom
  info_GPP = GOF(dobs.GPP_obs, dsim.GPP)
  info_ET = GOF(dobs.ET_obs, dsim.ET)

  if verbose
    indexes = [:KGE, :NSE, :R2, :RMSE, :bias, :bias_perc]
    info_GPP = info_GPP[indexes] |> round2
    info_ET = info_ET[indexes] |> round2
  end
  # _gpp = info_GPP[of_gof]
  _gpp = all(isnan.(dobs.GPP_obs)) ? 0.6 : info_GPP[of_gof]

  gof = [info_ET[of_gof], _gpp]
  weighted_mean(gof, [1.0, 1.0]) # goal
end


## gof of ET and GPP for ALL IGBP
function model_gof(df_out::DataFrame; all=true)
  lst = groupby(df_out, :IGBP)

  fun_et(d) = GOF(d.ET_obs, d.ET)
  fun_gpp(d) = GOF(d.GPP_obs, d.GPP)
  func = (; ET=fun_et, GPP=fun_gpp)

  map(f -> begin
      res = map_df_tuple(fun_et, lst)
      all && push!(res, (; IGBP="ALL", fun_et(df_out)...))
      DataFrame(res)
    end, func)
end
