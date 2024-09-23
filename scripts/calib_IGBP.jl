IGBPs = unique_sort(data.IGBPname)
IGBP = IGBPs[1]

df = data[data.IGBP.==IGBP, :]
printstyled("[IGBP=$IGBP] \n", bold=true, color=:green, underline=true)

# @time theta, goal, flag = m_calib(df; IGBPcode, maxn=2500);
@time _theta, goal, flag = model_calib(df, par0; IGBPcode=df.IGBPcode[1], maxn=2500);

model_goal(df, _theta; verbose=true)

r = PMLV2_sites(df; par=theta2par(_theta))
r = cbind(df[:, [:site, :date, :ET_obs, :GPP_obs]], r)

GOF(r.ET_obs, r.ET)
GOF(r.GPP_obs, r.GPP)
