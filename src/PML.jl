module PML

export PMLV2, PMLV2_sites, 
  photosynthesis, cal_Ei_Dijk2021, 
  T_adjust_Vm25, f_VPD_Zhang2019
export DataFrame, GOF
export file_FLUXNET_CRO

using DocStringExtensions
using Parameters, DataFrames
import HydroTools: cal_Uz, ET0_eq, Cp, atm, GOF, sceua

## global data
dir_proj = "$(@__DIR__)/.."
file_FLUXNET_CRO = "$dir_proj/data/CRO/FLUXNET_CRO" |> abspath


# lambda: [MJ kg-1]
W2mm(Ra; λ) = Ra * 86400 / 1e6 / λ

include("main_Ipaper.jl")
include("Params.jl")
include("Utilize/Utilize.jl")
include("ET_helper.jl")
# include("PET_equilibrium.jl")
# include("Ei_EvapIntercepted.jl")
# include("Ec_CanopyTrans.jl")
# include("Es_EvapSoil.jl")
include("water_constrain.jl")
include("photosynthesis.jl")
include("model_PMLV2.jl")

end # module PML
