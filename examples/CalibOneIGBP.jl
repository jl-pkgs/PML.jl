using PML, Ipaper, Test, RTableTools
df_out, df, par = deserialize(file_FLUXNET_CRO)
df.GPP_obs = df.GPPobs
df.ET_obs = df.ETobs
r = PMLV2(df; par)

# GOF(fval_soil, df.f_value_soil).MAE
# fval_soil - df.f_value_soil

# ## 模型参数率定
# @testset "CHECK PMLV2 RESULT" 
begin
  @test GOF(df.Eeq, r.Eeq).MAE <= 0.01
  @test GOF(df_out.Es_eq, r.Es_eq).MAE <= 0.002
  @test GOF(df_out.Ga, r.Ga).MAE <= 1e-10
  @test GOF(df_out.Gc, r.Gc_w).MAE <= 1e-3
  @test GOF(df_out.Ei, r.Ei).MAE <= 1e-8 # Ei passed Test
  @test GOF(df_out.Ec, r.Ec).MAE <= 0.002
  @test GOF(df_out.Es, r.Es).MAE <= 0.01
  @test GOF(df_out.Ecr, r.Ecr).MAE <= 0.002
  @test GOF(df_out.Eca, r.Eca).MAE <= 0.002

  @test GOF(df_out.ET_sim, r.ET).MAE <= 0.015
  @test GOF(df_out.GPP_sim, r.GPP).MAE <= 1E-8
end

begin
  using Plots
  gr(framestyle=:box, titlefontsize=12)
  
  var = :ET_sim
  plot(
    plot(df_out[:, :GPP_sim] - r[:, :GPP], label="GPP"),
    plot(df_out[:, :Gc] - r[:, :Gc_w], label="Gc"),
    plot(df_out[:, :Ga] - r[:, :Ga], label="Ga", ylims=(-1,1).*1e-6),
    
    plot(df.Eeq - r[:, :Eeq], label="Eeq"),
    plot(df_out[:, :Es_eq] - r[:, :Es_eq], label="Es_eq"),
    plot(df_out[:, :Eca] - r[:, :Eca], label="Eca"),
    plot(df_out[:, :Ecr] - r[:, :Ecr], label="Ecr"),
    plot(df_out[:, :Ec] - r[:, :Ec], label="Ec"),
    plot(df_out[:, :Es] - r[:, :Es], label="Es"),
    plot(df.f_value_soil - r.fval_soil, label="fval_soil"),

    plot(df_out[:, :Ei] - r[:, :Ei], label="Ei"),
    plot(df_out[:, :ET_sim] - r[:, :ET], label="ET"),
    size = (1000, 700)
  )
end

begin
  plot()
  plot!(r[:, :ET], label="Julia")
  plot!(df_out[:, :ET_sim], label="MATLAB")
end

begin
  inds=1:1000
  plot()
  plot!(r.fval_soil[inds], label="Julia")
  plot!(df.f_value_soil[inds], label="MATLAB")
end
