# include("leaf_BigLeaf.jl")
# include("leaf_TwoBigLeaf.jl")
# include("leaf_TwoLeaf.jl")
using Parameters

abstract type AbstractLeaf{FT<:AbstractFloat} end


"""
    struct BigLeaf{FT<:AbstractFloat} <: AbstractLeaf{FT}  

G_s = gs_sunlit * LAI_sunlit + gs_shaded * L_shaded
T = T(G_s)
"""
@with_kw mutable struct BigLeaf{FT} <: AbstractLeaf{FT}
  Lai::FT = 0.0 | "-"    # leaf area index
  gs::FT = 0.0  | "m s-1"    # stomatal conductance for h2o
  # Rs_c::FT = 0.0
  # Rs_s::FT = 0.0
  Rn_c::FT = 0.0
  Rn_s::FT = 0.0
  # H_s::FT = 0.0   # _s: soil
  # LE_s::FT = 0.0
  # H_c::FT = 0.0   # _c: canopy
  T::FT = 0.0 | "mm d-1"  # Transpiration
  GPP::FT = 0.0 | "gC m-2 d-1"
end

function radiative_transfer!(leaf::BigLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}
  (; Lai) = leaf
  τ = exp(-kA * Lai) # 辐射冠层截留的比例
  leaf.Rn_c = (1 - τ) * Rn
  leaf.Rn_s = τ * Rn
end


# 导度不进行积分，直接乘到冠层，陈镜明BEPS
"""
T = T(gs_sunlit) * LAI_sunlit + T(gs_shaded) * Lai_shaded
"""
@with_kw mutable struct TwoLeaf{FT} <: AbstractLeaf{FT}
  Lai_sunlit::FT = 0.0
  Lai_shaded::FT = 0.0
  T_sunlit::FT = 0.0
  T_shaded::FT = 0.0
  gs_sunlit::FT = 0.0 # stomatal conductance for h2o
  gs_shaded::FT = 0.0 # stomatal conductance for h2o
  GPP_sunlit::FT = 0.0
  GPP_shaded::FT = 0.0
end

# 导度积分到冠层，戴永久CLM
"""
T = T(Gs_sunlit) + T(Gs_shaded)
"""
@with_kw mutable struct TwoBigLeaf{FT} <: AbstractLeaf{FT}
  Lai_sunlit::FT = 0.0
  Lai_shaded::FT = 0.0
  gs_sunlit::FT = 0.0 # stomatal conductance for h2o
  gs_shaded::FT = 0.0 # stomatal conductance for h2o
  T_sunlit::FT = 0.0
  T_shaded::FT = 0.0
  GPP_sunlit::FT = 0.0
  GPP_shaded::FT = 0.0
end


@with_kw mutable struct Leaves{FT} <: AbstractLeaf{FT}
  nlyr::Int = 10

  Lai_sunlit::Vector{FT} = zeros(FT, nlyr)
  Lai_shaded::Vector{FT} = zeros(FT, nlyr)

  Tleaf_sunlit::Vector{FT} = zeros(FT, nlyr)
  Tleaf_shaded::Vector{FT} = zeros(FT, nlyr)

  gs_sunlit::Vector{FT} = zeros(FT, nlyr) # stomatal conductance for h2o
  gs_shaded::Vector{FT} = zeros(FT, nlyr) # stomatal conductance for h2o

  Rs::Vector{FT} = zeros(FT, nlyr + 1) # the last is ground
  Rl::Vector{FT} = zeros(FT, nlyr + 1)
  Rn::Vector{FT} = zeros(FT, nlyr + 1)

  LE_sunlit::Vector{FT} = zeros(FT, nlyr)
  LE_shaded::Vector{FT} = zeros(FT, nlyr)

  GPP_sunlit::FT = zeros(FT, nlyr)
  GPP_shaded::FT = zeros(FT, nlyr)
end


# BigLeaf{Float32}()
# Ec_CanopyTrans
# Es_EvapSoil
