using Ipaper, HydroTools, PML
using DataFrames, RTableTools
using Test, BenchmarkTools

function replace_miss!(df::AbstractDataFrame, miss=NaN)
  # colnames = names(df)
  # num_cols = [name for name in colnames if getDataType(df[!, name]) <: Number]
  for col in names(df)
    x = df[!, col]
    type = getDataType(x)

    if type <: AbstractFloat
      # replace_miss(x, miss)
      df[!, col] = drop_missing(x)
    end
  end
  df
end

f = "data/PMLv2_training_forcing_flux_v20200828 (80%)_102sp.csv"
data = fread(f)
replace_miss!(data)

## 处理观测值
data.GPP_obs = nanmean2.(data.GPP_NT, data.GPP_DT)
data.ET_obs = W2mm.(data.LE, data.Tavg)
data.Rn = cal_Rn.(data.Rs, data.Rl_in, data.Tavg, data.Albedo, data.Emiss)

data = @rename(data, Ca = CO2) # LAI_raw, LAI_whit, LAI_sgfitw
data.LAI = data.LAI_sgfitw
