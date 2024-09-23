using PML, RTableTools, HydroTools

df = fread("data/CRO/INPUT.csv")
df.Ca = df.CO2

df_out = fread("data/CRO/OUTPUT.csv")
param = fread("data/CRO/PARAM.csv")

inds = [1:4; 10:11; 5:9]
_names = names(param)[inds]
_values = Matrix(param)[inds]

par = Param_PMLV2(_values..., 0.5)
r = PMLV2_sites(df; par)
GOF(df_out.ETsim, r.ET)

GOF(df_out.ETsim, r.ET)
GOF(df_out.ETsim, r.ET)

# 蒸发部分，仍然不能完全匹配上，仍需要核对
df_out.ETsim - r.ET
df_out.GPPsim - r.GPP
df_out.ETobs - r.ET

begin
  using Plots
  gr(framestyle=:box)
  plot()
  # plot(df_out.ETobs, label="OBS")
  plot!(df_out.ETsim, label="MATLAB")
  plot!(r.ET, label="Julia")
end

begin
  plot(df_out.GPPsim, label="MATLAB")
  plot!(r.GPP, label="Julia")
end
