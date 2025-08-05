using SPAC, Test
import DataFrames: DataFrame

# include("test-Ipaper.jl")
include("test-stomatal_conductance.jl")
include("test-photosynthesis.jl")
include("test-PET.jl")

include("test-radiation.jl")
include("test-evapotranspiration.jl")


@testset "Model Params update!" begin
  FT = Float64
  model = Photosynthesis_Rong2018{FT}()

  params = Params(model)
  params |> DataFrame

  parnames = [:kQ, :VCmax25, :VPDmin]
  parvalues = [0.6, 10., 0.8]
  @time update!(model, parnames, parvalues; params)

  @test model.kQ == 0.6
  @test model.VCmax25 == 10.
  @test model.watercons.VPDmin == 0.8
end


FT = Float64
model = LandModel{FT}(;
  stomatal=Stomatal_Yu2004{FT}(),
  photosynthesis=Photosynthesis_Rong2018{FT}()
)
Params(model) |> DataFrame
