function model_goal(df, theta; IGBPcode=nothing, of_gof=:NSE, verbose=false)
  # IGBPcode !== nothing && (par.hc = hc_raw[IGBPcode])
  IGBPcode !== nothing && (theta[end] = hc_raw[IGBPcode]) # the last one is hc
  par = Param_PMLV2(theta...)

  dobs = df[!, [:GPP_obs, :ET_obs]]
  dsim = PMLV2_sites(df; par)

  ## 8-day, yearly, yearly_anom
  info_GPP = GOF(dobs.GPP_obs, dsim.GPP)
  info_ET = GOF(dobs.ET_obs, dsim.ET)

  if verbose
    indexes = [:KGE, :NSE, :R2, :RMSE, :bias, :bias_perc]
    info_GPP = info_GPP[indexes] |> round2
    info_ET = info_ET[indexes] |> round2

    @show info_GPP
    @show info_ET
  end

  goal = (info_ET[of_gof] + info_GPP[of_gof]) / 2
  goal
end


## 最后一步，参数率定模块
"""
    model_calib(df::AbstractDataFrame, par0::AbstractETParam; 
        IGBPcode=nothing, maxn=2500, of_gof=:NSE, kw...)
"""
function model_calib(df::AbstractDataFrame, par0::AbstractETParam; IGBPcode=nothing, maxn=2500, of_gof=:NSE, kw...)
  parRanges = get_bounds(par0)
  lower = parRanges[:, 1]
  upper = parRanges[:, 2]
  theta0 = collect(par0) # must be vector

  sceua(theta -> -model_goal(df, theta; IGBPcode, of_gof, kw...),
    theta0, lower, upper; maxn, kw...) # theta, goal, flag
end

export model_goal, model_calib
