using PenmanMonteithLeuning, Test
import HydroTools: atm

@testset "photosynthesis" begin
  # Prcp = 2.0 # mm
  # Rn = 50.0
  # U2 = 2.0
  Tavg = 20.0
  Rs = 200.0
  VPD = 2.0
  LAI = 2.0
  Ca = 380 # ppm

  par = par0
  GPP, Gc_w = photosynthesis(Tavg, Rs, VPD, LAI, atm, Ca; par)
  rs = 1 / Gc_w # s m-1, 阻力
  @test GPP ≈ 8.622226891575519
  @test Gc_w ≈ 0.0021840273263829023 # m s-1
  @test rs ≈ 457.86972897274114
end
# 1 / (0.46 * mol2m(Tavg))
