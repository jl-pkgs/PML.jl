
begin
  plot(d.LAI_whit, label="WHIT")
  plot!(d.LAI_glass, label="GLASS")
end

run_model(; site)
