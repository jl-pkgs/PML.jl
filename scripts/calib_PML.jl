
function calib_PML(data::AbstractDataFrame; of_gof=:KGE, maxn=2500)
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

# data.LAI .= data.LAI_raw
# data.LAI .= data.LAI_sgfitw
# data.LAI .= data.LAI_whit

# data = @rename(data, LAI = LAI_whit) # LAI_raw, LAI_whit, LAI_sgfitw
r = calib_PML(data; of_gof=:KGE, maxn=500)
r.gof

# out = fread("D:/GitHub/PML/PMLV2_Kong2019.m/OUTPUT/PMLv2_flux102_Cal_flux_v012.csv")
# out = @rename(out, ET = ETsim, GPP = GPPsim, ET_obs = ETobs, GPP_obs = GPPobs)
# model_gof(out)

## LAI_raw
#     NSE       R2        KGE
# ET  0.659226  0.681471  0.823243
# GPP 0.724708  0.730096  0.796629

## LAI_sgfitw
#     NSE       R2        KGE
# ET  0.659226  0.681471  0.823243
# GPP 0.724708  0.730096  0.796629

## LAI_whit
#     NSE       R2        KGE
# ET  0.613259  0.646577  0.803562
# GPP 0.676971  0.682734  0.78375
