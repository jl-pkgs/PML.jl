# include("leaf_BigLeaf.jl")
# include("leaf_TwoBigLeaf.jl")
# include("leaf_TwoLeaf.jl")


abstract type AbstractLeaf{FT<:AbstractFloat} end

"""
    struct BigLeaf{FT<:AbstractFloat} <: AbstractLeaf{FT}  
"""
@with_kw mutable struct BigLeaf{T} <: AbstractLeaf{FT}
  LAI::FT = 1.0
  Ec::FT = 1.0
  GPP::FT = 1.0
end

# 积分到冠层
@with_kw mutable struct TwoBigLeaf{T} <: AbstractLeaf{FT}
  LAI_sunlit::FT = 1.0
  LAI_shaded::FT = 1.0
  Ec_sunlit::FT = 1.0
  Ec_shaded::FT = 1.0
  GPP_sunlit::FT = 1.0
  GPP_shaded::FT = 1.0
end

# 不进行积分，直接乘到冠层
@with_kw mutable struct TwoLeaf{T} <: AbstractLeaf{FT}
  LAI_sunlit::FT = 1.0
  LAI_shaded::FT = 1.0
  Ec_sunlit::FT = 1.0
  Ec_shaded::FT = 1.0
  GPP_sunlit::FT = 1.0
  GPP_shaded::FT = 1.0
end
