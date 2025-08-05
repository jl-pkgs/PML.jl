export to_mat;


## DATATYPE CONVERSION ---------------------------------------------------------
function to_mat(res::output_PML{T}) where {T<:Real}
  names = fieldnames(output_PML)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data, names
end

function to_df(res::output_PML{T}) where {T<:Real}
  data, names = to_mat(res)
  DataFrame(data, names)
end
