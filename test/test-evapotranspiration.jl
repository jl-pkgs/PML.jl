using SPAC, Test


@testset "evapotranspiration" begin
  T = Float64

  air = AirLayer{T}()
  canopy_bare = BigLeaf{T}(; Lai = 0.0)
  canopy_veg = BigLeaf{T}(; Lai = 2.0)
  evap = Evapotranspiration_PML{T}()
  photo = Photosynthesis_Rong2018{T}()
  stomatal = Stomatal_Yu2004{T}()
  stomatal = Stomatal_Medlyn2011{T}()

  r_bare = evapotranspiration(air, canopy_bare, evap, photo, stomatal)
  r_veg = evapotranspiration(air, canopy_veg, evap, photo, stomatal)

  r_bare.ET_water == r_veg.ET_water
  @test r_bare[[:GPP, :Ec, :Ecr, :Eca, :Ei, :Pi]] == (; GPP=0.0, Ec=0.0, Ecr=0.0, Eca=0.0, Ei=0.0, Pi=0.0)
end
