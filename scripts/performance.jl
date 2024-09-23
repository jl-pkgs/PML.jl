of_gof = :KGE; maxn = 1000
@profview params = par_map(IGBP -> begin
    df = data[data.IGBP.==IGBP, vars]
    IGBPcode = df.IGBPcode[1]
    theta, _, _ = model_calib(df, par0; IGBPcode, of_gof, maxn, verbose=false)
    theta
  end, IGBPs)
