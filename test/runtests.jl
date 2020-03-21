using Test
using SignalAnalysis
using SignalAnalysis.Units

@testset "signals" begin

  x = signal(randn(8000), 1000)
  @test x isa AbstractArray
  @test length(x) == 8000
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 1
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000

  x = signal(randn((8000,2)), 1000)
  @test x isa AbstractArray
  @test size(x) == (8000, 2)
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 2
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000

  x1 = signal(x, 1000)
  @test x === x1

  x1 = signal(x, 8000)
  @test framerate(x) == 1000
  @test framerate(x1) == 8000

  x = @rate 1000 randn(8000)
  @test x isa AbstractArray
  @test length(x) == 8000
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 1
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000

  x1 = @samerateas x randn(8000)
  @test x1 isa AbstractArray
  @test nframes(x1) == 8000
  @test nchannels(x1) == 1
  @test framerate(x1) == 1000

  x1 = @samerateas 1//2 x randn(8000)
  @test x1 isa AbstractArray
  @test nframes(x1) == 8000
  @test nchannels(x1) == 1
  @test framerate(x1) == 500

  x1 = signal(randn(8000), 1000)
  x = analytic(x1)
  @test nframes(x1) == nframes(x)
  @test nchannels(x1) == nchannels(x)
  @test framerate(x1) == framerate(x)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x2 = analytic(x)
  @test isanalytic(x2)
  @test x === x2

  x1 = signal(randn((8000,2)), 1000)
  x = analytic(x1)
  @test nframes(x1) == nframes(x)
  @test nchannels(x1) == nchannels(x)
  @test framerate(x1) == framerate(x)
  @test isanalytic(x)
  @test real(x) ≈ x1
  x2 = analytic(x)
  @test isanalytic(x2)
  @test x === x2

  x = randn((8000,2))
  x1 = signal(x, 1000)
  @test x1 !== x
  @test samples(x1) === x
  @test samples(x) === x

  x = randn(8000)
  x1 = padded(x, (10, 5); delay=1)
  @test x1[-9] == 0
  @test x1[8005] == 0
  @test x1[1] == 0
  @test x1[2] == x[1]
  @test x1[8001] == x[8000]

  x1 = padded(x, (10, 5); delay=-2)
  @test x1[-9] == 0
  @test x1[8005] == 0
  @test x1[-1] == x[1]
  @test x1[7998] == x[8000]

  x = signal(ones(1000), 8kHz)
  x2 = slide(Float32, x, 250) do x1, blknum, firstframe
    sum(x1)*blknum
  end
  @test x2 == [250.0, 500.0, 750.0, 1000.0]
  x2 = slide(Tuple{Int,Float64}, x, 250) do x1, blknum, firstframe
    (blknum, sum(x1)*blknum)
  end
  @test x2 == [(1, 250.0), (2, 500.0), (3, 750.0), (4, 1000.0)]
  slide(x, 250) do x1, blknum, firstframe
    x1 .= blknum
  end
  @test (x[1], x[251], x[501], x[751]) == (1.0, 2.0, 3.0, 4.0)

  x = signal(randn(2000), 8kHz)
  @test toframe(0.2s, x) == 1601
  @test toframe([0.2s, 201ms], x) == [1601, 1609]
  @test toframe(0.2:0.01:0.3, x)[1:3] == [1601, 1681, 1761]

  x = signal(randn(8000), 1000)
  t = toframe(0:0.1:1, x)
  @test t == 1:100:1001

end

@testset "generate" begin

  x = cw(5000, 0.2, 44100)
  @test nframes(x) == 8820
  @test nchannels(x) == 1
  @test framerate(x) == 44100
  @test isanalytic(x)
  @test abs.(x) ≈ ones(8820)
  @test mean(ifrequency(x)) ≈ 5000

  x1 = cw(-5000, 0.2, 44100)
  @test samples(x) .* samples(x1) ≈ ones(8820)
  @test x .* x1 ≈ ones(8820)

  x1 = cw(5kHz, 200ms, 44.1kHz)
  @test samples(x1) ≈ samples(x)
  @test framerate(x1) == framerate(x)

  x1 = cw(5kHz, 200ms, 44.1kHz; window=(tukey, 0.05))
  @test nframes(x1) == 8820
  @test nchannels(x1) == 1
  @test framerate(x1) == 44100
  @test isanalytic(x1)
  @test energy(x1) < energy(x)
  @test abs(x1[1]) < abs(x1[2]) < 1
  @test abs(x1[end]) < abs(x1[end-1]) < 1

  x2 = cw(5kHz, 200ms, 44.1kHz; window=hamming)
  @test nframes(x2) == 8820
  @test nchannels(x2) == 1
  @test framerate(x2) == 44100
  @test isanalytic(x2)
  @test energy(x2) < energy(x1)
  @test abs(x2[1]) < abs(x2[2]) < 1
  @test abs(x2[end]) < abs(x2[end-1]) < 1

  x1 = cw(5kHz, 200ms, 44.1kHz; phase=90°)
  @test nframes(x1) == 8820
  @test nchannels(x1) == 1
  @test framerate(x1) == 44100
  @test isanalytic(x1)
  @test energy(x1) ≈ energy(x)
  @test samples(x1) .* -1im ≈ samples(x)

  x = chirp(5kHz, 7kHz, 100ms, 44.1kHz)
  @test nframes(x) == 4410
  @test nchannels(x) == 1
  @test framerate(x) == 44100
  @test isanalytic(x)
  @test abs.(x) ≈ ones(4410)

  x3 = ifrequency(x)
  @test x3[2] ≈ 5000
  @test x3[end-1] ≈ 7000

  x1 = chirp(7kHz, 5kHz, 100ms, 44.1kHz)
  @test nframes(x) == 4410
  @test nchannels(x) == 1
  @test framerate(x) == 44100
  @test isanalytic(x)
  @test abs.(x) ≈ ones(4410)

  x1 = ifrequency(x1)
  @test x1[3] ≈ 7000 atol=10
  @test x1[end-2] ≈ 5000 atol=10

  x1 = chirp(5kHz, 7kHz, 100ms, 44.1kHz; shape=:hyperbolic, window=(tukey,0.05))
  @test nframes(x1) == 4410
  @test nchannels(x1) == 1
  @test framerate(x1) == 44100
  @test isanalytic(x1)
  @test energy(x1) < energy(x)
  @test abs(x1[1]) < abs(x1[2]) < 1
  @test abs(x1[end]) < abs(x1[end-1]) < 1

  x2 = ifrequency(x1)
  @test x2[3] ≈ 5000 atol=10
  @test x2[2205] < x3[2205]
  @test x2[end-2] ≈ 7000 atol=10

  x2 = chirp(5kHz, 7kHz, 100ms, 44.1kHz; phase=45°, window=hamming)
  @test nframes(x2) == 4410
  @test nchannels(x2) == 1
  @test framerate(x2) == 44100
  @test isanalytic(x2)
  @test energy(x2) < energy(x1)
  @test abs(x2[1]) < abs(x2[2]) < 1
  @test abs(x2[end]) < abs(x2[end-1]) < 1

  x2 = chirp(5kHz, 7kHz, 100ms, 44.1kHz; phase=90°)
  @test nframes(x2) == 4410
  @test nchannels(x2) == 1
  @test framerate(x2) == 44100
  @test isanalytic(x2)
  @test energy(x2) ≈ energy(x)
  @test samples(x2) .* -1im ≈ samples(x)

  x1 = chirp(7kHz, 5kHz, 100ms, 44.1kHz; shape=:hyperbolic)
  @test nframes(x1) == 4410
  @test nchannels(x1) == 1
  @test framerate(x1) == 44100
  @test isanalytic(x1)

  x2 = ifrequency(x1)
  @test x2[2] ≈ 7000
  @test x2[end-1] ≈ 5000

end

@testset "basic" begin

  x1 = signal(randn(8000), 1000)
  x2 = signal(randn((8000,2)), 1000)

  @test energy(x1) isa Real
  @test energy(x1) ≈ sum(abs2, x1)/1000
  @test energy(samples(x1)) isa Real
  @test energy(samples(x1)) ≈ 1000*energy(x1)
  @test energy(samples(x1); fs=1000) ≈ energy(x1)

  @test size(energy(x2)) == (2,)
  @test energy(x2) ≈ [energy(x2[:,1]), energy(x2[:,2])]
  @test energy(samples(x2)) ≈ 1000*energy(x2)
  @test energy(samples(x2); fs=1000) ≈ energy(x2)

  x1 = signal(vcat(ones(8000), zeros(4000)), 1000)
  x2 = signal(ones((8000,2)), 1000)

  @test meantime(x1) ≈ 4 atol=0.01
  @test meantime(samples(x1); fs=1000) ≈ 4 atol=0.01
  @test meantime(x2) ≈ [4, 4] atol=0.01
  @test meantime(samples(x2); fs=1000) ≈ [4, 4] atol=0.01

  @test rmsduration(x1) ≈ 2.3094 atol=0.01
  @test rmsduration(samples(x1); fs=1000) ≈ 2.3094 atol=0.01
  @test rmsduration(x2) ≈ [2.3094, 2.3094] atol=0.01
  @test rmsduration(samples(x2); fs=1000) ≈ [2.3094, 2.3094] atol=0.01

  x1 = cw(5000, 0.2, 44100)
  x2 = [x1 x1]
  @test mean(ifrequency(x1)) ≈ 5000
  @test meanfrequency(x1) ≈ 5000 atol=100
  @test dropdims(mean(ifrequency(x2); dims=1); dims=1) ≈ [5000, 5000]
  @test meanfrequency(x2) ≈ [5000, 5000] atol=100

  x1 = chirp(5kHz, 7kHz, 100ms, 44.1kHz)
  x2 = [x1 chirp(7kHz, 5kHz, 100ms, 44.1kHz)]
  @test meanfrequency(x1) ≈ 6000 atol=100
  @test meanfrequency(x2) ≈ [6000, 6000] atol=100*√2
  @test rmsbandwidth(x1) ≈ 577.35 atol=100
  @test rmsbandwidth(x2) ≈ [577.35, 577.35]  atol=100

end

@testset "dsp" begin

  x = 1 .+ randn(100000)

  f = fir(127, 0, 5kHz; fs=44.1kHz)
  @test length(f) == 127
  y = filt(f, 1, x)
  @test energy(y) < energy(x)*2/3
  @test mean(y) ≈ mean(x) atol=0.1

  f = fir(127, 15kHz; fs=44.1kHz)
  @test length(f) == 127
  y = filt(f, 1, x)
  @test energy(y) < energy(x)/2
  @test mean(y) ≈ 0 atol=0.01

  f = fir(127, 5kHz, 10kHz; fs=44.1kHz)
  @test length(f) == 127
  y = filt(f, 1, x)
  @test energy(y) < energy(x)/2
  @test mean(y) ≈ 0 atol=0.01

  @test mean(removedc(x)) ≈ 0 atol=0.01
  @test mean(x) ≈ 1 atol=0.01
  removedc!(x)
  @test mean(x) ≈ 0 atol=0.01

  x = randn(100000) .* (1 .+ sin.(2π*50 .* (1:100000)./44100))
  x = signal(x, 44.1kHz)
  y = demon(x)
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10
  y = demon(samples(x); fs=framerate(x))
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10

  x = randn(100000,2) .* (1 .+ sin.(2π*50 .* (1:100000)./44100))
  x = signal(x, 44.1kHz)
  y = demon(x)
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2
  y = demon(samples(x); fs=framerate(x))
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2

end

@testset "rand" begin

  x = rand(RedGaussian(100000))
  @test length(x) == 100000
  @test mean(x) ≈ 0 atol=1e-1
  @test std(x) ≈ 1 atol=1e-1
  x1 = 0
  x2 = 0
  for j ∈ 1:100
    x1 += abs2(sum(x .* cis.(2π*(100+j) .* (1:100000)./10000)))
    x2 += abs2(sum(x .* cis.(2π*(1000+j) .* (1:100000)./10000)))
  end
  @test 1 < log10(x1/x2) < 3

  x = rand(PinkGaussian(100000))
  @test length(x) == 100000
  @test mean(x) ≈ 0 atol=1e-1
  @test std(x) ≈ 1 atol=1e-1
  x1 = 0
  x2 = 0
  for j ∈ 1:100
    x1 += abs2(sum(x .* cis.(2π*(100+j) .* (1:100000)./10000)))
    x2 += abs2(sum(x .* cis.(2π*(1000+j) .* (1:100000)./10000)))
  end
  @test 0 < log10(x1/x2) < 2

end
