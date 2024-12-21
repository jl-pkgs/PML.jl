using DataFrames

function replace_miss!(df::AbstractDataFrame)
  # colnames = names(df)
  # num_cols = [name for name in colnames if getDataType(df[!, name]) <: Number]
  for col in names(df)
    x = df[!, col]
    type = getDataType(x)
    if type <: AbstractFloat
      df[!, col] = drop_missing(x)
    end
  end
  df
end
