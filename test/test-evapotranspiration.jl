using SPAC, Test


@testset "evapotranspiration" begin
  FT = Float64

  air = AirLayer{FT}()
  canopy_bare = BigLeaf{FT}(; Lai = 0.0)
  canopy_veg = BigLeaf{FT}(; Lai = 2.0)
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()
  stomatal = Stomatal_Medlyn2011{FT}()

  r_bare = evapotranspiration(air, canopy_bare, evap, photo, stomatal)
  r_veg = evapotranspiration(air, canopy_veg, evap, photo, stomatal)

  r_bare.ET_water == r_veg.ET_water
  @test r_bare[[:GPP, :Ec, :Ecr, :Eca, :Ei, :Pi]] == zeros(6)

  model = LandModel{FT}(; evap, photo, stomatal)
  Params(model) |> DataFrame
end
