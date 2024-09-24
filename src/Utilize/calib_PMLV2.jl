# all sites
function calib_PMLV2(data::AbstractDataFrame; of_gof=:KGE, maxn=2500)
  vars = ["IGBPname", "IGBPcode", "site", "date", "GPP_obs", "ET_obs",
    "Prcp", "Tavg", "U2", "Rn", "Rs", "VPD", "LAI", "Pa", "Ca"]
  IGBPs = unique_sort(data.IGBPname)

  @time params = par_map(IGBP -> begin
      df = data[data.IGBP.==IGBP, vars]
      IGBPcode = df.IGBPcode[1]
      theta, _, _ = model_calib(df, par0; IGBPcode, of_gof, maxn, verbose=false)
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

export calib_PMLV2
