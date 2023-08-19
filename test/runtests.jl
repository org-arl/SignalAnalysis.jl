using Test, Statistics, LinearAlgebra, DSP, DSP.Util
using Plots
using SignalAnalysis
using SignalAnalysis.Units
using StableRNGs

rng = StableRNG(0)

@testset "signals" begin

  x = signal(randn(rng, 8000), 1000)
  @test x isa AbstractArray
  @test length(x) == 8000
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 1
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000
  ndx = rand(1:nframes(x))
  @test time(ndx, x) == domain(x)[ndx]

  x = signal(randn(rng, (8000,2)), 1000)
  @test x isa AbstractArray
  @test size(x) == (8000, 2)
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 2
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000
  ndx = rand(1:nframes(x))
  @test time(ndx, x) == domain(x)[ndx]

  x1 = signal(x, 1000)
  @test x === x1

  x1 = signal(x, 8000)
  @test framerate(x) == 1000
  @test framerate(x1) == 8000

  x = signal(randn(rng, 8000), 1000)
  @test x isa AbstractArray
  @test length(x) == 8000
  @test nframes(x) == 8000
  @test framerate(x) == 1000
  @test nchannels(x) == 1
  @test !isanalytic(x)
  @test domain(x) ≈ 0:1/1000:7999/1000

  x1 = signal(randn(rng, 8000), framerate(x))
  @test x1 isa AbstractArray
  @test nframes(x1) == 8000
  @test nchannels(x1) == 1
  @test framerate(x1) == 1000

  x1 = signal(randn(rng, 8000), framerate(x)/2)
  @test x1 isa AbstractArray
  @test nframes(x1) == 8000
  @test nchannels(x1) == 1
  @test framerate(x1) == 500

  x1 = signal(randn(rng, 8000), 1000)
  x = analytic(x1)
  @test nframes(x1) == nframes(x)
  @test nchannels(x1) == nchannels(x)
  @test framerate(x1) == framerate(x)
  @test isanalytic(x)
  @test √2*real(x) ≈ x1
  @test rms(x) ≈ rms(x1) atol=0.001
  x2 = analytic(x)
  @test isanalytic(x2)
  @test x === x2

  x1 = signal(randn(rng, (8000,2)), 1000)
  x = analytic(x1)
  @test nframes(x1) == nframes(x)
  @test nchannels(x1) == nchannels(x)
  @test framerate(x1) == framerate(x)
  @test isanalytic(x)
  @test √2*real(x) ≈ x1
  @test rms(x) ≈ rms(x1) atol=0.001
  x2 = analytic(x)
  @test isanalytic(x2)
  @test x === x2

  x = randn(rng, (8000,2))
  x1 = signal(x, 1000)
  @test x1 !== x
  @test samples(x1) === x
  @test samples(x) === x

  x = randn(rng, 8000)
  x1 = padded(x, (10, 5); delay=1)
  @test x1[-9] == 0
  @test x1[8005] == 0
  @test x1[1] == 0
  @test x1[2] == x[1]
  @test x1[8001] == x[8000]
  x2 = padded(x1, (1, 1); delay=0)
  @test x2[-10] == 0
  @test x2[8006] == 0
  @test x2[1] == 0
  @test x2[8001] == x1[8001] == x[8000]
  x2d = randn(rng, 8000, 2)
  x2d1 = padded(x2d, (10, 5); delay=1)
  @test all(x2d1[-9,:] .== 0)
  @test all(x2d1[8005,:] .== 0)
  @test all(x2d1[1,:] .== 0)
  @test x2d1[2,:] == x2d[1,:]
  @test x2d1[8001,:] == x2d[8000,:]
  x2d2 = padded(x2d1, (1, 1); delay=0)
  @test all(x2d2[-10,:] .== 0)
  @test all(x2d2[8006,:] .== 0)
  @test x2d2[2,:] == x2d1[2,:] == x2d[1,:]
  @test x2d2[8001,:] == x2d1[8001,:] == x2d[8000,:]

  x1 = padded(x, (10, 5); delay=-2)
  @test x1[-9] == 0
  @test x1[8005] == 0
  @test x1[-1] == x[1]
  @test x1[7998] == x[8000]
  x2d1 = padded(x2d, (10, 5); delay=-2)
  @test all(x2d1[-9,:] .== 0)
  @test all(x2d1[8005,:] .== 0)
  @test x2d1[-1,:] == x2d[1,:]
  @test x2d1[7998,:] == x2d[8000,:]

  x = signal(collect(1:10), 1.0)
  @test length(collect(partition(x, 5))) == 2
  @test length(collect(partition(x, 4))) == 3
  @test length(collect(partition(x, 5; flush=false))) == 2
  @test length(collect(partition(x, 4; flush=false))) == 2
  @test length(collect(partition(x, 5; step=2))) == 5
  @test length(collect(partition(x, 4; step=2))) == 5
  @test length(collect(partition(x, 5; step=2, flush=false))) == 3
  @test length(collect(partition(x, 4; step=2, flush=false))) == 4

  pad_x = padded(x, (1, 1))
  @test iterate(partition(pad_x, 5)) == ([0,1,2,3,4], 5)
  @test iterate(partition(pad_x, 5), 1) == ([1,2,3,4,5], 6)
  @test length(collect(partition(pad_x, 5))) == 3
  @test length(collect(partition(pad_x, 4))) == 3
  @test length(collect(partition(pad_x, 5; flush=false))) == 2
  @test length(collect(partition(pad_x, 4; flush=false))) == 3
  @test length(collect(partition(pad_x, 5; step=2))) == 6
  @test length(collect(partition(pad_x, 4; step=2))) == 6
  @test length(collect(partition(pad_x, 5; step=2, flush=false))) == 4
  @test length(collect(partition(pad_x, 4; step=2, flush=false))) == 5

  x = signal(ones(1000), 8kHz)
  x2 = map(enumerate(partition(x, 250))) do (blknum, x1)
    @test size(x1) == (250,)
    @test framerate(x) == framerate(x1)
    sum(x1)*blknum
  end
  @test x2 == [250.0, 500.0, 750.0, 1000.0]
  x2 = map(enumerate(partition(x, 250))) do (blknum, x1)
    @test size(x1) == (250,)
    (blknum, sum(x1)*blknum)
  end
  @test x2 == [(1, 250.0), (2, 500.0), (3, 750.0), (4, 1000.0)]
  x2 = hcat(map(partition(x, 250)) do x1
    @test size(x1) == (250,)
    [sum(x1), prod(x1)]
  end...)'
  @test x2 == [250.0 1.0; 250.0 1.0; 250.0 1.0; 250.0 1.0;]
  for (blknum, x1) ∈ enumerate(partition(x, 250))
    x1 .= blknum
  end
  @test (x[1], x[251], x[501], x[751]) == (1.0, 2.0, 3.0, 4.0)

  x = signal(ones(1000, 2), 8kHz)
  x2 = map(enumerate(partition(x, 250))) do (blknum, x1)
    @test size(x1) == (250, 2)
    sum(x1)*blknum
  end
  @test x2 == [2*250.0, 2*500.0, 2*750.0, 2*1000.0]
  for x1 ∈ partition(x, 250)
    @test size(x1) == (250, 2)
  end

  x = signal(randn(rng, 2000), 8kHz)
  @test toframe(0.2s, x) == 1601
  @test toframe([0.2s, 201ms], x) == [1601, 1609]
  @test toframe((0.2, 0.201), x) == 1601:1609
  @test toframe(0.2:0.201s, x) == 1601:1609
  @test toframe(0.2:0.01:0.3, x)[1:3] == [1601, 1681, 1761]
  @test toframe(0.2s, framerate(x)) == 1601
  @test toframe([0.2s, 201ms], framerate(x)) == [1601, 1609]
  @test toframe((0.2, 0.201), framerate(x)) == 1601:1609
  @test toframe(0.2:0.201s, framerate(x)) == 1601:1609
  @test toframe(0.2:0.01:0.3, framerate(x))[1:3] == [1601, 1681, 1761]
  @test x[0.2:0.201s] == x[1601:1609]
  @test x[0:0.201s] == x[1:1609]

  x = signal(randn(rng, 8000), 1000)
  t = toframe(0:0.1:1, x)
  @test t == 1:100:1001
  t = toframe(0:0.1:1, framerate(x))
  @test t == 1:100:1001

  xlen = 8000
  x = randn(rng, xlen)
  fs = 1000
  s0 = signal(x, fs)
  @test vec(s0) == s0
  @test framerate(vec(s0)) == fs
  @test size(vec(s0)) == (xlen,)
  s1 = signal(reshape(x, :, 1), fs)
  @test vec(s1) == s0
  @test framerate(vec(s1)) == fs
  @test size(vec(s1)) == (xlen,)
  s2 = signal(randn(rng, xlen, 2), fs)
  @test_throws ArgumentError vec(s2)
  s3 = signal(randn(rng, xlen, 1, 1), fs)
  @test_throws ArgumentError vec(s3)

  x = signal(randn(rng, 100), 1000)
  x1 = samerateas(x)(randn(rng, 200))
  @test framerate(x1) == framerate(x)
  @test nframes(x1) == 200
  x1 = samerateas(x, randn(rng, 200))
  @test framerate(x1) == framerate(x)
  @test nframes(x1) == 200
  x1 = samerateas(randn(rng, 100), randn(rng, 200))
  @test framerate(x1) == 1.0
  @test nframes(x1) == 200

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

  x1 = signal(randn(rng, 8000), 1000)
  x2 = signal(randn(rng, (8000,2)), 1000)

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

  x = 1 .+ randn(rng, 100000)

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

  x = randn(rng, 100000) .* (1 .+ sin.(2π*50 .* (1:100000)./44100))
  x = signal(x, 44.1kHz)
  y = demon(x)
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10
  y = demon(x; method=:mean)
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10
  y = demon(x; method=:fir)
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10
  y = demon(samples(x); fs=framerate(x))
  @test length(y) == 400
  @test meanfrequency(y) ≈ 50 atol=10

  x = randn(rng, 100000,2) .* (1 .+ sin.(2π*50 .* (1:100000)./44100))
  x = signal(x, 44.1kHz)
  y = demon(x)
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2
  y = demon(x; method=:mean)
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2
  y = demon(x; method=:fir)
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2
  y = demon(samples(x); fs=framerate(x))
  @test size(y) == (400, 2)
  @test meanfrequency(y) ≈ [50, 50] atol=10*√2

  x = signal(2*round.(Int,rand(rng, 1000)).-1 + 1im*(2*round.(Int,rand(rng, 1000)).-1), 100)/√2
  y = upconvert(x, 10, 100)
  @test framerate(y) == 1000
  @test rms(y[100:end-100])*sqrt(10) ≈ 1.0 atol=0.01
  @test meanfrequency(y) ≈ 100.0 atol=5.0
  z = downconvert(y, 10, 100)[12:end-11]
  @test framerate(z) == 100
  @test amp2db(rms(x-z)/rms(z)) < -40.0

  x = signal(2*round.(Int,rand(rng, 1000,2)).-1 + 1im*(2*round.(Int,rand(rng, 1000,2)).-1), 100)/√2
  y = upconvert(x, 10, 100)
  @test nchannels(y) == 2
  @test framerate(y) == 1000
  @test rms(y[100:end-100,1])*sqrt(10) ≈ 1.0 atol=0.01
  @test rms(y[100:end-100,2])*sqrt(10) ≈ 1.0 atol=0.01
  @test meanfrequency(y) ≈ [100.0, 100.0] atol=5.0
  z = downconvert(y, 10, 100)[12:end-11,:]
  @test nchannels(z) == 2
  @test framerate(z) == 100
  @test amp2db(rms(x-z)/rms(z)) < -40.0

  x = x[:,1:1]
  y = upconvert(x, 10, 100)
  @test nchannels(y) == 1
  @test framerate(y) == 1000
  @test y isa AbstractMatrix
  z = downconvert(y, 10, 100)[12:end-11,:]
  @test nchannels(z) == 1
  @test framerate(z) == 100
  @test z isa AbstractMatrix

  x = signal(rand(rng, [-1,1], 1000), 1000.0)
  @test upconvert(x, 8, 2000.0) == upconvert(complex.(x), 8, 2000.0)
  x = signal(rand(rng, [-1,1], 1000, 2), 1000.0)
  @test upconvert(x, 8, 2000.0) == upconvert(complex.(x), 8, 2000.0)

  x = rcosfir(0.25, 10)
  @test length(x) == 221
  @test x[1:10:end] ≈ vcat(zeros(11), 0.326598866403003, zeros(11))
  @test sum(x.^2) ≈ 1.0

  x = rrcosfir(0.25, 10)
  @test length(x) == 221
  y = conv(x, x)[1:10:end]
  @test y ≈ vcat(zeros(22), 1.0, zeros(22)) atol=0.01
  @test sum(x.^2) ≈ 1.0

  for j ∈ 2:10
    x = mseq(j)
    @test x isa Array{Float64}
    @test length(x) == 2^j-1
    @test all(abs.(x) .== 1)
    y = circconv(x)
    @test all(y .== circconv(x, x))
    @test y[1] ≈ length(x)
    @test all(y[2:end] .≈ -1)
    x = gmseq(j)
    @test x isa Array{ComplexF64}
    @test length(x) == 2^j-1
    @test all(abs.(x) .== 1)
    y = circconv(x)
    @test y[1] ≈ length(x)
    @test rms(y[2:end]) ≈ 0 atol=1e-10
  end

  @test mseq((1,3)) == mseq(3)
  @test all(gmseq(3, 0) .== 1.0)

  x = cw(10kHz, 0.1s, 80kHz)
  @test abs(goertzel(x, 10kHz)) ≈ length(x) atol=1e-6
  @test abs(goertzel(x, 9kHz)) ≈ 0 atol=1e-6
  @test abs(goertzel(x, 11kHz)) ≈ 0 atol=1e-6
  @test abs.(goertzel(x, 10kHz, 512))[1:end-1] ≈ repeat([512.0], 15) atol=1e-6
  @test abs.(goertzel(x, 9375, 512))[1:end-1] ≈ repeat([0.0], 15) atol=1e-6
  x = hcat(cw(10kHz, 0.1s, 80kHz), cw(11kHz, 0.1s, 80kHz))
  @test abs.(goertzel(x, 10kHz)) ≈ [size(x,1), 0] atol=1e-6
  @test abs.(goertzel(x, 9kHz)) ≈ [0, 0] atol=1e-6
  @test abs.(goertzel(x, 11kHz)) ≈ [0, size(x,1)] atol=1e-6

  x = cw(1kHz, 1s, 80kHz)
  x2 = x .+ 0.25*rand(rng, size(x)...) .+ 0.25im*rand(rng, size(x)...)
  e2 = √mean(abs2.(x[800:end] .- x2[800:end]))
  x3 = pll(x2)
  @test meanfrequency(x3[800:end]) ≈ meanfrequency(x[800:end]) atol=0.1
  e3 = √mean(abs2.(x[800:end] .- x3[800:end]))
  @test e2/e3 > 3

  x = cw(7kHz, 1s, 44.1kHz)
  f = fir(127, 5kHz, 10kHz; fs=44.1kHz)
  x1 = sfilt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = filt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = sfilt(f, 1, x)
  @test framerate(x1) == framerate(x)
  x1 = filt(f, 1, x)
  @test framerate(x1) == framerate(x)
  x1 = sfiltfilt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = filtfilt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = sresample(x, 3//2)
  @test framerate(x1) == 3 * framerate(x) / 2
  x1 = resample(x, 3//2)
  @test framerate(x1) == 3 * framerate(x) / 2

  x = signal(randn(rng, 100), 10kHz)
  x1 = signal(vcat(zeros(1000), x/2, zeros(1000)), 10kHz)
  x2 = mfilter(x, x1)
  @test !isanalytic(x2)
  @test argmax(abs.(x2)) == 1001

  x = analytic(signal(randn(rng, 100), 10kHz))
  x1 = signal(vcat(zeros(1000), x/2, zeros(1000)), 10kHz)
  x2 = mfilter(x, x1)
  @test isanalytic(x2)
  @test argmax(abs.(x2)) == 1001

  x = analytic(signal(randn(rng, 100), 10kHz))
  x1 = signal(randn(rng, 1000), 10kHz)
  x2 = mfilter(x, x1)
  @test isanalytic(x2)

  x = signal(randn(rng, 100), 10kHz)
  x1 = analytic(signal(randn(rng, 1000), 10kHz))
  x2 = mfilter(x, x1)
  @test isanalytic(x2)

  x = signal(randn(rng, 100), 10kHz)
  x1 = signal(randn(rng, 1000), 11kHz)
  @test_throws ArgumentError mfilter(x, x1)

  x = signal(randn(rng, 100), 10kHz)
  x1 = randn(rng, 1000)
  x2 = mfilter(x, x1)
  @test !isanalytic(x2)

  x = randn(rng, 100)
  x1 = signal(randn(rng, 1000), 10kHz)
  x2 = mfilter(x, x1)
  @test !isanalytic(x2)

  onesideds = [true, false]
  nfft = 256
  windows = [nothing, rect, tukey(nfft, 0.5), hanning, hamming]
  noverlaps = [0, 0, (15*nfft)÷16, (5*nfft)÷6, (5*nfft)÷6]
  for onesided in onesideds
    if onesided === true
      x = randn(rng, 96000)
    else
      x = randn(rng, 96000) + im .* randn(rng, 96000)
    end
    for (window, noverlap) in zip(windows, noverlaps)
      xstft = stft(x, nfft, noverlap; window=window, onesided=onesided)
      outputtype = isreal(x) ? Real : Complex
      x̂ = istft(outputtype, xstft; nfft=nfft, noverlap=noverlap, window=window)
      outputtype === Complex && (@test x̂ == istft(xstft; nfft=nfft, noverlap=noverlap, window=window))
      @test rms(x[nfft:length(x̂)-nfft]-x̂[nfft:end-nfft]) / rms(x[nfft:length(x̂)-nfft]) ≈ 0. atol=0.001
    end
  end

  fs = 9600
  hpf = fir(9, 1000; fs=fs)
  x = filtfilt(hpf, randn(rng, 96000))
  y = 0.1 .* real(chirp(500, 1000, 1.0, fs))
  x[9600+1:2*9600] .+= y
  nfft = 256
  noverlap = 0
  x̃ = whiten(x; nfft=nfft, noverlap=noverlap, window=nothing)
  ỹ = x̃[9600+1:2*9600]
  @test cor(y, ỹ) > 0.3
  x̃ = whiten(x; nfft=nfft, noverlap=noverlap, window=nothing, γ=0.)
  @test rms(x[1:length(x̃)] - x̃) / rms(x[1:length(x̃)]) ≈ 0. atol=0.0001

  x = chirp(1000, 5000, 0.1, 40960; window=(tukey, 0.05))
  x4 = resample(x, 4)
  y4 = samerateas(x4, zeros(32768))
  y4[128:127+length(x4)] = real(x4)          # time 0.000775s, index 32.75
  y4[254:253+length(x4)] += -0.8 * real(x4)  # time 0.001544s, index 64.25
  y4[513:512+length(x4)] += 0.6 * real(x4)   # time 0.003125s, index 129.0
  y = resample(y4, 1//4)
  y .+= 0.1 * randn(rng, length(y))
  t, a = findsignal(x, y, 3; coarse=false)
  @test t ≈ [0.000775, 0.001545, 0.003124] atol=2e-6
  @test real(a) / real(a[1]) ≈ [1.0, -0.8, 0.6] atol=1e-2
  t, a = findsignal(x, y, 3; coarse=true)
  @test t ≈ [0.000775, 0.001545, 0.003124] atol=1e-5
  @test real(a) / real(a[1]) ≈ [1.0, -0.8, 0.6] atol=1e-2

  @test delay!([1,2,3,4], -1) == [2,3,4,0]
  @test delay!([1,2,3,4], 1) == [0,1,2,3]

  x = signal(vcat(zeros(128), [1.0, 2.0, 3.0, 2.0, 1.0], zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=0.05
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(128), ComplexF64[1.0, 2.0, 3.0, 2.0, 1.0], zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=0.05
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(ComplexF32, 128), ComplexF32[1.0, 2.0, 3.0, 2.0, 1.0], zeros(ComplexF32, 128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=0.05
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

end

@testset "rand" begin

  x = rand(rng, RedGaussian(100000))
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

  x = rand(rng, PinkGaussian(100000))
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

  types = [Float16, Float32, Float64]
  for t ∈ types
    @test eltype(rand(rng, RedGaussian(1000, one(t)))) == t
    @test eltype(rand(rng, PinkGaussian(1000, one(t)))) == t
  end

end

@testset "plots" begin

  x = chirp(5kHz, 7kHz, 500ms, 44.1kHz)

  # plots may perhaps be tested with VisualRegressionTests.jl
  # for now we just test that the recipes don't die
  p = plot(x)
  @test p isa Plots.Plot
  p = plot(real(x))
  @test p isa Plots.Plot
  p = plot([x -x])
  @test p isa Plots.Plot
  p = plot(real([x -x]))
  @test p isa Plots.Plot
  p = plot(signal(randn(rng, 61040, 2), 1000))
  @test p isa Plots.Plot
  p = psd(x)
  @test p isa Plots.Plot
  p = psd(real(x))
  @test p isa Plots.Plot
  p = psd(real(x); xscale=:log10)
  @test p isa Plots.Plot
  p = psd([x -x])
  @test p isa Plots.Plot
  p = psd(real([x -x]))
  @test p isa Plots.Plot
  p = psd(real([x -x]); xscale=:log10)
  @test p isa Plots.Plot
  p = psd(samples(x); fs=framerate(x))
  @test p isa Plots.Plot
  p = psd(samples([x -x]); fs=framerate(x))
  @test p isa Plots.Plot
  p = specgram(x)
  @test p isa Plots.Plot
  p = specgram(x; pooling=nothing)
  @test p isa Plots.Plot
  p = specgram(samples(x); fs=framerate(x))
  @test p isa Plots.Plot
  f = fir(127, 15kHz; fs=44.1kHz)
  p = plotfreqresp(f)
  @test p isa Plots.Plot
  p = plotfreqresp(f; xscale=:log10)
  @test p isa Plots.Plot
  p = plotfreqresp(f, [1])
  @test p isa Plots.Plot

end

@testset "array" begin

  c = 1500.0
  θ = range(0.0, π; length=181)
  sd1 = steering(0.0:1.0:5.0, c, θ)
  @test size(sd1) == (6, 181)
  @test sd1[:,1] ≈ (2.5:-1.0:-2.5)/c atol=1e-6
  @test sd1[:,91] ≈ zeros(6) atol=1e-6
  @test sd1[:,181] ≈ (-2.5:1.0:2.5)/c atol=1e-6
  sd2 = steering(-2.0:1.0:3.0, c, θ)
  @test sd1 == sd2
  sd2 = steering(vcat((0.0:1.0:5.0)', zeros(1,6)), c, θ)
  @test sd1 == sd2
  sd2 = steering(vcat((0.0:1.0:5.0)', zeros(2,6)), c, θ)
  @test sd1 == sd2

  θ = range(0.0, 2π; length=9)
  rxpos = [-1.0 0.0 1.0 0.0; 0.0 1.0 0.0 -1.0]
  sd1 = steering(rxpos, c, θ)
  @test size(sd1) == (4, 9)
  @test sd1[:,1] ≈ [1.0, 0.0, -1.0, 0.0]/c atol=1e-6
  @test sd1[:,3] ≈ [0.0, -1.0, 0.0, 1.0]/c atol=1e-6
  @test sd1[:,5] ≈ -[1.0, 0.0, -1.0, 0.0]/c atol=1e-6
  @test sd1[:,7] ≈ -[0.0, -1.0, 0.0, 1.0]/c atol=1e-6
  @test sd1[1,4] ≈ sd1[2,4] atol=1e-6
  @test sd1[3,4] ≈ sd1[4,4] atol=1e-6
  @test sd1[3,4]-sd1[2,4] ≈ √2/c atol=1e-6

  θ = deg2rad.(1:360)
  fc = 500.0
  fs = 44100.0
  sd1 = steering(rxpos, c, θ)
  @test size(sd1) == (4, 360)
  x = cw(fc, 1.0, fs)
  x̄ = padded(x, 0; delay=round(Int, √2/c*fs))
  x4 = [x x x̄ x̄]
  bfo = beamform(x4, fc, 4096, sd1)
  @test argmax(bfo) == 135
  bfo = beamform(samples(x4), fc, 4096, sd1; fs=framerate(x4))
  @test argmax(bfo) == 135
  bfo = beamform(x4, fc, 4096, sd1; method=Bartlett())
  @test argmax(bfo) == 135
  bfo = beamform(x4 .+ 0.001*randn(rng, size(x4)), fc, 4096, sd1; method=Capon())
  @test argmax(bfo) == 135
  bfo = beamform(x4, fc, 4096, sd1; method=Capon(0.1))
  @test argmax(bfo) == 135
  bfo = beamform(x4, fc, 4096, sd1; method=Music())
  @test argmax(bfo) == 135
  bfo = beamform(x4, fc, 4096, sd1; method=Music(1))
  @test argmax(bfo) == 135
  y4 = goertzel(x4, fc, 4096)
  bfo = beamform(y4, fc, sd1)
  @test argmax(bfo) == 135
  bfo = beamform(x4, sd1)
  @test nchannels(bfo) == 360
  e = energy(bfo)
  @test argmax(e) == 135
  @test maximum(e) ≈ 16.0 atol=0.1
  @test meanfrequency(bfo[:,135]) ≈ fc atol=10.0

end

@testset "tfa" begin

  x = real.(chirp(5kHz, 10kHz, 1s, 44.1kHz))
  y = tfd(x, Spectrogram())
  @test y isa SignalAnalysis.TFD
  @test time(y) ≈ 0.0029024943310657597:0.005804988662131519:0.9955555555555555
  @test freq(y) ≈ range(0.0, 22050.0; length=129)
  @test power(y) isa AbstractMatrix
  @test size(power(y)) == (129, 172)
  @test freq(y)[argmax(power(y)[:,86])] ≈ 7500.0 atol=100.0
  y = tfd(x, Spectrogram(nfft=512, noverlap=256, window=hamming))
  @test y isa SignalAnalysis.TFD
  @test time(y) ≈ 0.005804988662131519:0.005804988662131519:0.9926530612244898
  @test freq(y) ≈ range(0.0, 22050.0; length=257)
  @test power(y) isa AbstractMatrix
  @test size(power(y)) == (257, 171)
  @test freq(y)[argmax(power(y)[:,86])] ≈ 7500.0 atol=10.0

  x = real.(chirp(5kHz, 10kHz, 0.01s, 44.1kHz))
  y = tfd(x, Wigner())
  @test y isa SignalAnalysis.TFD
  @test time(y) ≈ range(0.0, 0.00998; length=441)
  @test freq(y) ≈ range(0.0, 22050.0, length=442)
  @test power(y) isa AbstractMatrix
  @test size(power(y)) == (442, 441)
  @test freq(y)[argmax(power(y)[:,220])] ≈ 7500.0 atol=1.0
  y = tfd(x, Wigner(nfft=512, smooth=10, method=:CM1980, window=hamming))
  @test y isa SignalAnalysis.TFD
  @test time(y) ≈ range(0.0, 0.00998; length=441)
  @test freq(y) ≈ range(0.0, 22050.0, length=257)
  @test power(y) isa AbstractMatrix
  @test size(power(y)) == (257, 441)
  @test freq(y)[argmax(power(y)[:,220])] ≈ 7500.0 atol=100.0

end
