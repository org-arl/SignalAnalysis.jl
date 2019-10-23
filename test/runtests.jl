using Test
using SignalAnalysis
using AxisArrays

@testset "Signal" begin

    x1 = randn(8000)
    x = sampledsignal(x1, 8000)
    @test x isa AxisArray
    @test ndims(x) == 1
    @test length(x) == 8000
    @test samplingrate(x) == 8000.0
    @test !isanalytic(x)
    @test x1 == x
    @test time(x) == 0:1/8000:7999/8000

    x1 = randn((8000,2))
    x = sampledsignal(x1, 8000, starttime=1.2)
    @test x isa AxisArray
    @test ndims(x) == 2
    @test size(x) == (8000,2)
    @test samplingrate(x) == 8000.0
    @test !isanalytic(x)
    @test x1 == x
    @test time(x) == 1.2 .+ (0:1/8000:7999/8000)

    x1 = randn(8000)
    x = sampledsignal(x1, 8000, analytic=true)
    @test x isa AxisArray
    @test ndims(x) == 1
    @test length(x) == 8000
    @test samplingrate(x) == 8000.0
    @test isanalytic(x)
    @test real(x) ≈ x1

    x = sampledsignal(x)
    @test x isa AxisArray
    @test ndims(x) == 1
    @test length(x) == 8000
    @test samplingrate(x) == 8000.0
    @test isanalytic(x)
    @test real(x) ≈ x1

    x = sampledsignal(x, analytic=false)
    @test x isa AxisArray
    @test ndims(x) == 1
    @test length(x) == 8000
    @test samplingrate(x) == 8000.0
    @test !isanalytic(x)
    @test x ≈ x1

    x1 = randn((8000,2))
    x = sampledsignal(x1, 8000, analytic=true)
    @test x isa AxisArray
    @test ndims(x) == 2
    @test size(x) == (8000,2)
    @test samplingrate(x) == 8000.0
    @test isanalytic(x)
    @test real(x) ≈ x1

    x = sampledsignal(x, analytic=false)
    @test x isa AxisArray
    @test ndims(x) == 2
    @test size(x) == (8000,2)
    @test samplingrate(x) == 8000.0
    @test !isanalytic(x)
    @test x ≈ x1

end

@testset "Basic" begin

    x1 = randn(8000)
    x = sampledsignal(x1, 8000)
    @test energy(x) isa Number
    @test energy(x) == energy(x1)
    @test energy(x) ≈ sum(abs2, x1)

    x1 = randn((8000,2))
    x = sampledsignal(x1, 8000)
    @test size(energy(x)) == (1,2)
    @test energy(x) == energy(x1)
    @test energy(x) ≈ sum(abs2, x1, dims=1)

end
