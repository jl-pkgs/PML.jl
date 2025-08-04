"""
    PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca; par=param0, frame=3)

# Notes
一个站点的计算。注意，不同植被类型，参数不同。

# Arguments
- `frame`: in 8-days
"""
function PMLV2(Prcp::V, Tavg::V, Rs::V, Rn::V,
  VPD::V, U2::V, LAI::V,
  Pa::V,
  Ca::V,
  PC::Union{T,V}=T(1.0);
  # par::LandModel, 
  frame=3,
  res::Union{Nothing,output_PML}=nothing) where {T<:Real,V<:AbstractVector{T}}

  n = length(Prcp)
  fields = fieldnames(interm_PML)[2:end-2] # skip ET, fval_soil and Es
  r = interm_PML{T}()
  res === nothing && (res = output_PML{T}(; n))

  # isvec_Ca = isa(Ca, AbstractVector)
  isvec_PC = isa(PC, AbstractVector)
  @inbounds for t = 1:n
    # _Ca = isvec_Ca ? Ca[t] : Ca
    _PC = isvec_PC ? PC[t] : PC

    PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t],
      Ca[t], _PC; par, r)
    res[t, fields] = r
  end

  res.fval_soil = movmean2(Prcp, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.fval_soil, 0, 1)

  res.Es .= res.fval_soil .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end


"""
# Arguments
- `kw`: named keyword arguments
  + `r`: `interm_PML`
"""
function PMLV2(d::AbstractDataFrame; par::LandModel, kw...)

  (; Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca) = d
  PC = "PC" ∈ names(d) ? d.PC : 1.0

  PMLV2(Prcp, Tavg, Rs, Rn,
    VPD, U2, LAI,
    Pa, Ca, PC; par, kw...) |> to_df
end


# 相同植被类型多个站点一起的计算
function PMLV2_sites(df::AbstractDataFrame; par::LandModel, kw...)
  "site" ∉ names(df) && return PMLV2(df; par, kw...)

  sites = df.site
  grps = unique(sites)

  # ! 这里有优化空间，但程序会写的很复杂
  res = []
  for grp in grps
    inds = sites .== grp
    d = df[inds, :]
    r = PMLV2(d; par, kw...)
    push!(res, r)
  end
  vcat(res...)
end
