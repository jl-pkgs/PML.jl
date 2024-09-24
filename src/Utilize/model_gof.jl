function map_df_tuple(fun::Function, lst::GroupedDataFrame{DataFrame}, args...; kw...)
  n = length(lst)
  _keys = keys(lst)
  map(i -> begin
      d = lst[i]
      key = NamedTuple(_keys[i])
      r = fun(d, args...; kw...)
      (; key..., r...)
    end, 1:n)
end

## 统计每种植被类型的模拟效果
function model_gof(df_out::DataFrame; all=true)
  fun_et(d) = GOF(d.ET_obs, d.ET)
  fun_gpp(d) = GOF(d.GPP_obs, d.GPP)

  lst = groupby(df_out, :IGBP)

  res = map_df_tuple(fun_et, lst)
  all && push!(res, (; IGBP="ALL", fun_et(df_out)...))
  gof_et = DataFrame(res)

  res = map_df_tuple(fun_gpp, lst)
  all && push!(res, (; IGBP="ALL", fun_gpp(df_out)...))
  gof_gpp = DataFrame(res)
  (; ET=gof_et, GPP=gof_gpp)
end

export model_gof, map_df_tuple
