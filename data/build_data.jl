using PML

param = fread("data/CRO/PARAM.csv")
df_out = fread("data/CRO/OUTPUT.csv")
df_out.Ec = df_out.Ecr + df_out.Eca

df = fread("data/CRO/INPUT.csv")
df.Ca = df.CO2

inds = [1:4; 10:11; 5:9]
_names = names(param)[inds]
_values = Matrix(param)[inds]
par = Param_PMLV2(_values..., 0.5)

l = (; df_out, df, par)
serialize("data/CRO/FLUXNET_CRO", l)
