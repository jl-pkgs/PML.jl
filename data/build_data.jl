using PML


function read_data()
  df_out = fread("data/CRO/OUTPUT.csv")
  df_out.Ec = df_out.Ecr + df_out.Eca

  df = fread("data/CRO/INPUT.csv")
  df.Ca = df.CO2

  df_out, df
end

param = fread("data/CRO/PARAM.csv")
inds = [1:4; 10:11; 5:9]
_names = Symbol.(names(param)[inds])
_values = Matrix(param)[inds]
_par = NamedTuple(_names .=> _values)
# par = Param_PMLV2(_values..., 0.5)

l = (; df_out, df, _par)
serialize("data/CRO/FLUXNET_CRO", l)

## all sites
site = "US-Twt"
inds = findall(df.site .== site)
d_out = df_out[inds, :]

_par = (α=0.03265625, η=0.069296875, g1=9.552734375, Am_25=17.671875, VPDmin=1.21515625, VPDmax=3.5, D0=0.6541015625, kQ=0.10114375, kA=0.89921875, S_sls=0.01015625, fER0=0.152734375, hc=0.5)

df_out, df = read_data()

l = (; df_out, df, _par)
serialize("data/CRO/FLUXNET_CRO_US-Twt", l)
