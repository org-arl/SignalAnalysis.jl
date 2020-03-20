using Test
using SignalAnalysis

@testset "Signal" begin

  x = signal(randn(8000), 8000)
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/8000:7999/8000

  x = signal(randn((8000,2)), 8000)
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/8000:7999/8000

  x1 = signal(randn(8000), 8000)
  x = analytic(x1)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x = analytic(x)
  @test isanalytic(x)
  @test real(x) ≈ x1

  x = signal(randn((8000,2)), 8000)
  x = analytic(x1)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x = analytic(x)
  @test isanalytic(x)
  @test real(x) ≈ x1

end

@testset "Basic" begin

  x = signal(randn(8000), 8000)
  @test energy(x) isa Real
  @test energy(x) ≈ sum(abs2, x)/8000

  x = randn((8000,2))
  @test size(energy(x)) == (2,)
  @test energy(x) ≈ [energy(x[:,1]), energy(x[:,2])]

end
