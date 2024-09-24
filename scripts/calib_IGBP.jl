IGBPs = unique_sort(data.IGBPname)
IGBP = IGBPs[1]

df = data[data.IGBP.==IGBP, :]
printstyled("[IGBP=$IGBP] \n", bold=true, color=:green, underline=true)

# @time theta, goal, flag = m_calib(df; IGBPcode, maxn=2500);
r = PMLV2_sites(df; par=theta2par(_theta))
r = cbind(df[:, [:site, :date, :ET_obs, :GPP_obs]], r)

GOF(r.ET_obs, r.ET)
GOF(r.GPP_obs, r.GPP)
