function map_df_tuple(fun::Function, lst::GroupedDataFrame{DataFrame}, args...; kw...)
  n = length(lst)
  _keys = keys(lst)
  map(i -> begin
      d = lst[i]
      key = NamedTuple(_keys[i])
      r = fun(d, args...; kw...)
      (; key..., r...)
    end, 1:n)
end


## 统计每种植被类型的模拟效果
function model_gof(df_out::DataFrame; all=true)
  fun_et(d) = GOF(d.ET_obs, d.ET)
  fun_gpp(d) = GOF(d.GPP_obs, d.GPP)

  lst = groupby(df_out, :IGBP)

  res = map_df_tuple(fun_et, lst)
  all && push!(res, (; IGBP="ALL", fun_et(df_out)...))
  gof_et = DataFrame(res)

  res = map_df_tuple(fun_gpp, lst)
  all && push!(res, (; IGBP="ALL", fun_gpp(df_out)...))
  gof_gpp = DataFrame(res)
  (; ET=gof_et, GPP=gof_gpp)
end


function model_goal(df, theta; IGBPcode=nothing, of_gof=:NSE, verbose=false)
  # IGBPcode !== nothing && (par.hc = hc_raw[IGBPcode])
  IGBPcode !== nothing && (theta[end] = hc_raw[IGBPcode]) # the last one is hc
  par = Param_PMLV2(theta...)

  dobs = df[!, [:GPP_obs, :ET_obs]]
  # dsim = PMLV2(df; par) # the exact MATLAB version
  dsim = PMLV2_sites(df; par) 

  ## 8-day, yearly, yearly_anom
  info_GPP = GOF(dobs.GPP_obs, dsim.GPP)
  info_ET = GOF(dobs.ET_obs, dsim.ET)

  if verbose
    indexes = [:KGE, :NSE, :R2, :RMSE, :bias, :bias_perc]
    info_GPP = info_GPP[indexes] |> round2
    info_ET = info_ET[indexes] |> round2
  end
  gof = [info_ET[of_gof], info_GPP[of_gof]]
  goal = weighted_mean(gof, [1, 0.5])
  goal
end


export model_goal
export model_gof, map_df_tuple
