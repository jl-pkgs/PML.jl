using PML, Test, Ipaper, DataFrames

df_out, df, par = deserialize(file_FLUXNET_CRO_USTwt)
r = PMLV2_sites(df; par)

@testset "PMLV2 scalar" begin
  Prcp = 2.0 # mm
  Tavg = 20.0
  Rs = 200.0
  Rn = 50.0
  VPD = 2.0
  U2 = 2.0
  LAI = 2.0
  @test_nowarn PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI;)
end

@testset "CHECK PMLV2 RESULT" begin
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


@testset "model_calib" begin
  df.GPP_obs = df.GPPobs
  df.ET_obs = df.ETobs
  @time _theta, goal, flag = model_calib(df, par0; IGBPcode=df.IGBPcode[1], maxn=2500)
  goal = model_goal(df, _theta; verbose=true)
  @test goal > 0.55 # mean(KGE_GPP, KGE_ET)
end

@testset "model_gof" begin
  d = DataFrame(; IGBP=["CRO", "CRO", "CRO"],
    ET_obs=[1.0, 2, 3], ET=[0.5, 0.6, 0.7],
    GPP_obs=[1.0, 2, 3], GPP=[1, 1.5, 2.0])
  @test size(model_gof(d).ET, 1) == 2
  @test size(model_gof(d; all=false).ET, 1) == 1
end
