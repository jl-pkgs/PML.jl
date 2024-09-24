module PML

export PMLV2_sites
export PMLV2, T_adjust_Vm25, f_VPD_Zhang2019

using DocStringExtensions
using Parameters, DataFrames
import HydroTools: cal_rho_a, cal_Uz, ET0_eq, 
  Cp, atm, 
  GOF, sceua

# lambda: [MJ kg-1]
W2mm(Ra; λ) = Ra * 86400 / 1e6 / λ

include("main_Ipaper.jl")
include("Params.jl")
include("Utilize/Utilize.jl")
include("ET_helper.jl")
# include("PET_equilibrium.jl")
include("Ei_EvapIntercepted.jl")
# include("Ec_CanopyTrans.jl")
# include("Es_EvapSoil.jl")
include("water_constrain.jl")
include("photosynthesis.jl")
include("model_PMLV2.jl")

end # module PML
