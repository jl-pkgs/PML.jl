module SPAC

export PMLV2, PMLV2_sites,
  photosynthesis, cal_Ei_Dijk2021,
  T_adjust_Vm25, f_VPD_Zhang2019

export file_FLUXNET_CRO, file_FLUXNET_CRO_USTwt
export DataFrame, GOF
export fread, fwrite, melt_list


using Parameters, FieldMetadata
using DataFrames
using Statistics
using RTableTools
using DocStringExtensions


# import HydroTools: GOF, sceua
import Ipaper: par_map


## global data
dir_proj = "$(@__DIR__)/.."
file_FLUXNET_CRO = "$dir_proj/data/CRO/FLUXNET_CRO" |> abspath
file_FLUXNET_CRO_USTwt = "$dir_proj/data/CRO/FLUXNET_CRO_US-Twt" |> abspath

include("ModelParam.jl")

include("modules/modules.jl")
# include("utilize.jl")

# include("ModelCalib.jl")
# include("PMLV2.jl")

# include("Ei_EvapIntercepted.jl")
# include("Ec_CanopyTrans.jl")
# include("Es_EvapSoil.jl")
# include("model_PMLV2.jl")

end # module SPACmodels
