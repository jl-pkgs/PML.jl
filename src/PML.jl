module PML

export PMLV2, PMLV2_sites, 
  photosynthesis, cal_Ei_Dijk2021, 
  T_adjust_Vm25, f_VPD_Zhang2019
export DataFrame, GOF
export file_FLUXNET_CRO, file_FLUXNET_CRO_USTwt

using Parameters
using FieldMetadata
using DataFrames
import HydroTools: cal_Uz, Cp, atm, GOF, sceua
using DocStringExtensions

## global data
dir_proj = "$(@__DIR__)/.."
file_FLUXNET_CRO = "$dir_proj/data/CRO/FLUXNET_CRO" |> abspath
file_FLUXNET_CRO_USTwt = "$dir_proj/data/CRO/FLUXNET_CRO_US-Twt" |> abspath


include("main_Ipaper.jl")
include("Params.jl")
include("Optim/ModelCalib.jl")
include("ET_helper.jl")
# include("PET_equilibrium.jl")
# include("Ei_EvapIntercepted.jl")
# include("Ec_CanopyTrans.jl")
# include("Es_EvapSoil.jl")
include("water_constrain.jl")
include("photosynthesis.jl")
include("model_PMLV2.jl")

end # module PML
