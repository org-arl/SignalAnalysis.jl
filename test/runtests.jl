using Test
using SignalAnalysis
using AxisArrays

@testset "Signal" begin

  samplingrate(8000)
  @test samplingrate() == 8000

  x = randn(8000)
  @test !isanalytic(x)
  @test time(x) ≈ 0:1/8000:7999/8000

  x = randn((8000,2))
  @test !isanalytic(x)
  @test time(x) ≈ 0:1/8000:7999/8000

  x1 = randn(8000)
  x = analytic(x1)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x = analytic(x)
  @test isanalytic(x)
  @test real(x) ≈ x1

  x = randn((8000,2))
  x = analytic(x1)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x = analytic(x)
  @test isanalytic(x)
  @test real(x) ≈ x1

end

@testset "Basic" begin

  samplingrate(8000)

  x = randn(8000)
  @test energy(x) isa Real
  @test energy(x) ≈ sum(abs2, x)/8000

  x = randn((8000,2))
  @test size(energy(x)) == (2,)
  @test energy(x) ≈ [energy(x[:,1]), energy(x[:,2])]

end
