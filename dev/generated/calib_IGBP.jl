using Plots
gr(framestyle=:box)
x = range(0, stop=6π, length=1000)
y1 = sin.(x)
y2 = cos.(x)
plot(x, [y1, y2])

using PML, Ipaper

df_out, df, par = deserialize(file_FLUXNET_CRO)
df.GPP_obs = df.GPPobs
df.ET_obs = df.ETobs

theta, goal, flag = model_calib(df, par0)
df_out = PMLV2_sites(df; par=theta2par(theta))

gof = [
  (; var="ET", GOF(df.ET_obs, df_out.ET)...),
  (; var="GPP", GOF(df.GPP_obs, df_out.GPP)...)] |> DataFrame
DataFrame(gof)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
