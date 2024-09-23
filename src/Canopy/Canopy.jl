# include("leaf_BigLeaf.jl")
# include("leaf_TwoBigLeaf.jl")
# include("leaf_TwoLeaf.jl")
using Parameters

abstract type AbstractLeaf{FT<:AbstractFloat} end

"""
    struct BigLeaf{FT<:AbstractFloat} <: AbstractLeaf{FT}  
"""
@with_kw mutable struct BigLeaf{FT} <: AbstractLeaf{FT}
  L::FT = 0.0
  gs::FT = 0.0 # stomatal conductance for h2o
  T::FT = 0.0
  GPP::FT = 0.0
end

# 导度积分到冠层，戴永久CLM
@with_kw mutable struct TwoBigLeaf{T} <: AbstractLeaf{FT}
  L_sunlit::FT = 0.0
  L_shaded::FT = 0.0
  gs_sunlit::FT = 0.0 # stomatal conductance for h2o
  gs_shaded::FT = 0.0 # stomatal conductance for h2o
  T_sunlit::FT = 0.0
  T_shaded::FT = 0.0
  GPP_sunlit::FT = 0.0
  GPP_shaded::FT = 0.0
end

# 导度不进行积分，直接乘到冠层，陈镜明BEPS
@with_kw mutable struct TwoLeaf{T} <: AbstractLeaf{FT}
  L_sunlit::FT = 0.0
  L_shaded::FT = 0.0
  T_sunlit::FT = 0.0
  T_shaded::FT = 0.0
  gs_sunlit::FT = 0.0 # stomatal conductance for h2o
  gs_shaded::FT = 0.0 # stomatal conductance for h2o
  GPP_sunlit::FT = 0.0
  GPP_shaded::FT = 0.0
end

# BigLeaf{Float32}()
