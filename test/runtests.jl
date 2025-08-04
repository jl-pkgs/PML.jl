using PenmanMonteithLeuning, Test
import DataFrames: DataFrame

include("test-Ipaper.jl")
include("test-stomatal_conductance.jl")

# include("test-PMLV2.jl")
# include("test-photosynthesis.jl")

@testset "ModelParams update!" begin
  FT = Float64
  model = Photosynthesis_Rong2018{FT}()

  params = ModelParams(model)
  params |> DataFrame

  parnames = [:kQ, :VCmax25, :VPDmin]
  parvalues = [0.6, 10., 0.8]
  @time update!(model, parnames, parvalues; params)

  @test model.kQ == 0.6
  @test model.VCmax25 == 10.
  @test model.watercons.VPDmin == 0.8
end
