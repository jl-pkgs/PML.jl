using RTableTools, Ipaper, PML
using Plots
include("../examples/main_pkgs.jl")

f = "Z:/Researches/ET_ModelDev/data-raw/Forcing_PMLV2_China_8day_2003-2022_flux37_v20250105_V2.csv"
df = fread(f) |> replace_miss!
sites = df.name |> unique_sort 
sites = setdiff(sites, ["元江"])

# site = "那曲"
site = "哀牢山"
d = df[df.name.==site, :]

parNames = [
  :α, :η, :g1, :Am_25, :VPDmin, :VPDmax, :D0, :kQ, :kA, :S_sls, :fER0#, :d_pc # :hc
]

import Base: NamedTuple
NamedTuple(names::AbstractVector, values::AbstractVector) = 
  NamedTuple{tuple(names...)}(values)

function run_model(; site="",
  use_PC=true, type_lai="whit", 
  of_gof=:NSE, kw...)

  d = df[df.name.==site, :]
  (; date, GPP, ET, prcp, LAI_whit, LAI_glass, Rnl, Rns, Rs, Tavg, U2, VPD, Ca, Pa, PC) = d

  Rn = Rns + Rnl
  forcing = (; date, GPP_obs=GPP, ET_obs=ET,
    Prcp=prcp, Tavg, Rs, Rn, VPD, U2,
    # LAI=LAI_glass, 
    # LAI=LAI_whit,
    Pa, Ca) |> DataFrame
  use_PC && (forcing.PC = PC)
  type_lai == "whit" && (forcing.LAI = LAI_whit)
  type_lai == "glass" && (forcing.LAI = LAI_glass)

  ## 模型参数率定
  theta, goal, flag = ModelCalib(forcing, par0, parNames; of_gof)
  par = theta2par(theta, parNames)

  out = PMLV2_sites(forcing; par)
  cbind(out; GPP_obs=GPP, ET_obs=ET)

  gof = [
    (; site, kw..., var="ET", GOF(ET, out.ET)...),
    (; site, kw..., var="GPP", GOF(GPP, out.GPP)...)] |> DataFrame
  (; out=out[:, Cols(:GPP_obs, :ET_obs, 1:end)], 
    par = NamedTuple(parNames, theta), gof)
end

# run_model(; site, use_PC=true, type_lai="whit", of_gof=:NSE)#.gof
# run_model(; site, use_PC=false, type_lai="whit", of_gof=:NSE).gof
function get_prefix(use_PC, type_lai)
  prefix = use_PC ? "WithPC" : "NonPC"
  if use_PC && :d_pc ∉ parNames
    prefix = "ConstPC"
  end
  "LAI_$type_lai" * ",$prefix"
end

function process(; use_PC=false, type_lai="glass")
  prefix = get_prefix(use_PC, type_lai)
  N = length(sites)
  res = Vector{Any}(undef, N)
  @par for i in 1:N
    site = sites[i]
    printstyled("[$i] $site\n", color=:green, bold=true)
    res[i] = run_model(; site, use_PC, type_lai, of_gof=:NSE)
  end

  df_gof = vcat(map(x -> x.gof, res)...)
  df_out = vcat(map(x -> x.out, res)...)
  
  mat_par = cat(map(x -> collect(x.par), res)..., dims=2) |> transpose |> collect
  df_par = cbind(DataFrame(mat_par, parNames); site = sites)[:, Cols(:site, 1:end)]

  fwrite(df_gof, "./OUTPUT/PMLV2China_flux37_$(prefix)_gof.csv")
  fwrite(df_par, "./OUTPUT/PMLV2China_flux37_$(prefix)_par.csv")
  fwrite(df_out, "./OUTPUT/PMLV2China_flux37_$(prefix)_OUTPUT.csv")

  GOF(df_out.GPP_obs, df_out.GPP)
  GOF(df_out.ET_obs, df_out.ET)
end

# process(; type_lai="whit", use_PC=false)
process(; type_lai="whit", use_PC=true)
process(; type_lai="glass", use_PC=false)
process(; type_lai="glass", use_PC=true)

# LAI_whit & Non_PC
NSE = 0.4901658678314159, R2 = 0.5713317663570666, KGE = 0.7435630375012554

# LAI_whit & With_PC
NSE = 0.5080601316516402, R2 = 0.5889340168907004, KGE = 0.7519358746627769

# LAI_glass & Non_PC
NSE = 0.4935561078799632, R2 = 0.5797932490279245, KGE = 0.7532742010794355

# LAI_glass & With_PC
NSE = 0.5343449820736267, R2 = 0.6058902430946861, KGE = 0.7663923308482895
