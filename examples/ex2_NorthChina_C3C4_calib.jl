using PML, Ipaper, Test, RTableTools

include("main_pkgs.jl")
data = fread("./Forcing_C3C4_sp2.csv")
replace_miss!(data)

df = data[data.site.=="固城", :]
# df = data[data.site.=="禹城", :]
# df.ET_obs = df.ET_obs ./ 0.9

par = par0
r = PMLV2(df; par)

theta, goal, flag = ModelCalib(df, par0)
df_out = PMLV2_sites(df; par=theta2par(theta))
# df_out[1:10, :]
sites = unique(df.site)


begin
  RES = []
  for site in sites
    inds = df.site .== site
    d_obs = df[inds, :]
    d_sim = df_out[inds, :]

    gof = [
      (; var="ET", GOF(d_obs.ET_obs, d_sim.ET)...),
      (; var="GPP", GOF(d_obs.GPP_obs, d_sim.GPP)...)] |> DataFrame
    gof = gof[:, [:var, :NSE, :KGE, :R2, :bias_perc, :n_valid]]
    push!(RES, gof)
  end
  RES
end

## PML-ET
# 固城: 0.46
# 禹城: 0.63
