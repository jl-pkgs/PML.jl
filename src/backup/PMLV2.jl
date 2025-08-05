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
