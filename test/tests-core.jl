function test_signals()

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

  x = signal(t -> sin(2π * 100 * t), 1, 3, 8000)
  @test x isa AbstractVector{Float64}
  @test nframes(x) == 16000
  @test framerate(x) == 8000
  @test x ≈ sin.(2π * 100 .* (1 .+ (0:15999) ./ 8000))
  @test signal(t -> sin(2π * 100 * t), 1𝓈, 3𝓈, 8kHz) == x
  x = signal(cispi, 2, 100)
  @test x isa AbstractVector{ComplexF64}
  @test nframes(x) == 200
  @test framerate(x) == 100
  @test x ≈ cispi.((0:199) ./ 100)
  @test signal(cispi, 2000ms, 0.1kHz) == x
  @test signal(cispi, 0, 2, 100) == x
  @test nframes(signal(sin, 100ms, 44.1kHz)) == 4410
  @test nframes(signal(sin, 1, 1, 100)) == 0
  @test @inferred(signal(sin, 0, 1, 100)) isa AbstractVector{Float64}
  @test @inferred(signal(sin, 1, 100)) isa AbstractVector{Float64}

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

  x = signal(collect(1.0:10.0), 1000)
  x1 = padded(x, (2, 3))
  @test axes(x1, 1) == 1:15
  @test x1[3] == x[1]
  @test x1[1] == 0
  @test framerate(x1) == 1000
  @test domain(x1) ≈ (0:14) ./ 1000
  @test samples(x1) == [0; 0; 1:10; 0; 0; 0]
  @test axes(padded(collect(1.0:10.0), (2, 3)), 1) == -1:13
  xr = padded(signal(randn(rng, 64), 8000), (2, 3))
  @test framerate(downconvert(xr, 1, 1000)) == 8000
  @test framerate(upconvert(downconvert(xr, 1, 1000), 1, 1000)) == 8000
  @test length(vcat(xr, xr)) == 138
  @test length(analytic(xr)) == 69
  @test ifft(fft(xr)) ≈ samples(xr)
  x1 = padded(x)
  @test axes(x1) == axes(x)
  @test length(x1) == 10
  @test x1[-100] == 0
  @test x1[1000] == 0
  @test x1[3] == x[3]
  @test samples(x1[-1:2]) == [0, 0, 1, 2]
  @test collect(view(x1, -1:2)) == [0, 0, 1, 2]
  @test framerate(x1) == 1000
  @test samples(x1) == samples(x)
  x1 = padded(x; delay=2, fill=-1)
  @test x1[1] == -1
  @test x1[3] == x[1]
  @test x1[12] == x[10]
  @test samples(x1) == [-1; -1; 1:8]
  x1 = padded(collect(1:3))
  @test x1[0] == 0 && x1[2] == 2
  x2d = padded(signal([1.0 2.0; 3.0 4.0], 1000))
  @test x2d[0, 2] == 0
  @test x2d[2, 2] == 4
  @test_throws BoundsError x2d[0, 3]
  @test_throws BoundsError padded(collect(1.0:4.0))[2, 5]
  @test padded(collect(1.0:4.0))[2, 1] == 2

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
  @test x[200:201ms] == x[1601:1609]
  @test x[(0.2s, 201ms)] == x[1601:1609]
  @test_throws ArgumentError x[(0.2, 0.201)]
  x[(0.2s, 0.201s)] = zeros(9)
  @test all(iszero, x[1601:1609])

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

  xc = signal(randn(rng, ComplexF64, 64), 1000)
  X = @inferred fft(xc)
  @test X isa SignalAnalysis.SampledSpectrum
  @test samples(X) isa Vector{ComplexF64}
  @test samples(X) ≈ fft(samples(xc))
  @test framerate(X) == 1000
  @test nframes(X) == 64
  @test @inferred(domain(X)) ≈ (0:63) .* (1000 / 64)
  @test typeof(domain(X)) == typeof(domain(xc))
  @test fftfreq(X) ≈ fftfreq(64, 1000)
  @test fftfreq(X) isa Vector{Float64}
  @test eltype(fftfreq(xc)) == eltype(freq(1:2, X)) == Float64
  @test eltype(domain(xc)) == Float64 && time(2, xc) isa Float64
  @test startswith(sprint(show, MIME("text/plain"), X), "SampledSpectrum of 64-point signal @ 1000.0 Hz, bins 0:63")
  @test occursin("one-sided", sprint(show, MIME("text/plain"), rfft(real(xc))))
  @test occursin("padded", sprint(show, MIME("text/plain"), padded(X)))
  y = @inferred ifft(X)
  @test y isa SignalAnalysis.SampledSignal
  @test framerate(y) == 1000
  @test samples(y) ≈ samples(xc)
  @test bfft(X) ≈ 64 .* samples(xc)
  @test framerate(bfft(X)) == 1000
  @test ifft(xc) ≈ ifft(samples(xc))
  @test bfft(xc) ≈ bfft(samples(xc))
  @test framerate(ifft(log.(abs.(X)))) == 1000

  for fs ∈ (8000, 11025, 44100, 48000, 96000), n ∈ (63, 100, 1000, 4097)
    x = signal(randn(rng, n), fs)
    @test framerate(ifft(fft(x))) == framerate(x)
    @test framerate(irfft(rfft(x))) == framerate(x)
  end

  B = X[11:20]
  @test B isa SignalAnalysis.SampledSpectrum
  @test domain(B) ≈ domain(X)[11:20]
  @test freq(1, B) ≈ domain(X)[11]
  @test freq(0, B) ≈ domain(X)[10]
  @test freq(1:3, B) ≈ domain(X)[11:13]
  @test domain(X[:]) == domain(X)
  @test domain(@view X[11:20]) ≈ domain(B)
  @test domain(X[1:4:end]) ≈ domain(X)[1:4:end]
  @test domain(X[5:4:end]) ≈ domain(X)[5:4:end]
  @test domain(X[11:1:20]) ≈ domain(B)
  @test_throws ArgumentError X[2:4:end]
  @test_throws ArgumentError X[1:3:end]
  @test_throws ArgumentError X[end:-1:1]
  @test X[[1, 3]] isa Vector{ComplexF64}
  @test X[abs.(samples(X)) .> 0] isa Vector{ComplexF64}
  @test X[3] == samples(X)[3]
  @test samples(ifft(X[1:4:end])) ≈ vec(sum(reshape(samples(xc), 16, 4); dims=2))

  @test_throws ArgumentError ifft(B)
  @test_throws ArgumentError irfft(B)
  @test X[125Hz:250Hz] isa SignalAnalysis.SampledSpectrum
  @test domain(X[125Hz:250Hz]) ≈ 125:15.625:250
  @test domain(X[125:250Hz]) ≈ domain(X[125Hz:250Hz])
  @test domain(X[(125Hz, 250Hz)]) ≈ domain(X[125Hz:250Hz])
  @test_throws ArgumentError X[-100Hz:100Hz]
  @test_throws ArgumentError X[(125, 250)]
  @test tobin(125Hz, X) == 9
  @test tobin(125, X) == 9
  @test tobin([125Hz, 250Hz], X) == [9, 17]
  @test tobin((125Hz, 250Hz), X) == 9:17
  @test tobin(125Hz, B) == -1
  @test tobin(-125Hz, X) == -7
  @test tobin(125Hz, 1kHz, 64) == 9
  @test tobin((125Hz, 250Hz), 1000, 64) == 9:17
  @test @inferred(tobin(125Hz, X)) isa Int
  Y = copy(X)
  Y[125Hz:250Hz] = zeros(9)
  @test all(iszero, samples(Y)[9:17])
  @test samples(Y)[18] == samples(X)[18]

  S = fftshift(X)
  @test samples(S) ≈ fftshift(samples(X))
  @test domain(S) ≈ (-32:31) .* (1000 / 64)
  @test fftfreq(S) ≈ domain(S)
  @test domain(S[-125Hz:125Hz]) ≈ (-8:8) .* (1000 / 64)
  @test_throws ArgumentError ifft(S)
  @test_throws ArgumentError fftshift(S)
  @test_throws ArgumentError ifftshift(X)
  @test_throws ArgumentError fftshift(B)
  @test_throws ArgumentError fftshift(rfft(real(xc)))
  @test ifftshift(S) ≈ X
  @test ifft(ifftshift(S)) ≈ samples(xc)
  x5 = signal(randn(rng, ComplexF64, 5), 1000)
  @test domain(fftshift(fft(x5))) ≈ (-2:2) .* 200
  @test ifft(ifftshift(fftshift(fft(x5)))) ≈ samples(x5)

  @test vcat(X[1:10], X[11:64]) ≈ X
  @test vcat(X[1:10], X[11:64]) isa SignalAnalysis.SampledSpectrum
  @test ifft(vcat(X[1:10], X[11:64])) ≈ samples(xc)
  @test vcat(X, ComplexF64[]) isa Vector{ComplexF64}
  @test vcat(ComplexF64[], X) isa Vector{ComplexF64}
  @test [zeros(4); S; zeros(4)] isa Vector{ComplexF64}
  @test_throws ArgumentError vcat(X[1:10], X[12:64])
  @test_throws ArgumentError vcat(X, X[1:1])
  @test_throws ArgumentError vcat(X[1:10], fft(signal(randn(rng, ComplexF64, 64), 2000))[11:64])
  @test hcat(X, X) isa SignalAnalysis.SampledSpectrum
  @test domain(hcat(B, B)) ≈ domain(B)
  @test hcat(X, samples(X)) isa Matrix{ComplexF64}
  @test_throws ArgumentError hcat(X, B)

  P = padded(B)
  @test P isa SignalAnalysis.SampledSpectrum
  @test axes(P) == axes(B)
  @test P[1:10] ≈ B
  @test P[0] == 0
  @test P[11] == 0
  @test P[1 + 64] == B[1]
  @test P[1 - 64] == B[1]
  @test domain(P[-5:15]) ≈ domain(X)[5:25]
  @test padded(P) === P
  @test_throws ArgumentError padded(P, 1)
  @test_throws ArgumentError fftshift(P)
  @test_throws ArgumentError padded(X; fill=1)
  @test_throws ArgumentError padded(B, 1; delay=1)
  Z = zeros(ComplexF64, 64)
  Z[11:20] .= samples(B)
  @test ifft(P) ≈ ifft(Z)
  @test framerate(ifft(P)) == 1000
  @test ifft(padded(S)) ≈ samples(xc)
  @test ifft(padded(S[-125Hz:125Hz])) ≈ ifft(samples(X) .* (abs.(fftfreq(64, 1000)) .≤ 125))
  @test domain(padded(X)[-125Hz:125Hz]) ≈ (-8:8) .* (1000 / 64)

  Q = padded(B, (2, 3))
  @test axes(Q, 1) == 1:15
  @test Q[3:12] ≈ B
  @test Q[1] == 0
  @test domain(Q) ≈ domain(X)[9:23]
  @test samples(Q) isa Vector{ComplexF64}
  @test vcat(Q, X[24:64]) isa SignalAnalysis.SampledSpectrum
  @test_throws ArgumentError padded(B, (30, 30))
  @test ifft(padded(X[2:64], (1, 0))) ≈ ifft([0; samples(X)[2:64]])

  xr = signal(randn(rng, 64), 1000)
  @test fft(xr) ≈ fft(samples(xr))
  R = @inferred rfft(xr)
  @test R ≈ rfft(samples(xr))
  @test domain(R) ≈ (0:32) .* (1000 / 64)
  @test fftfreq(R) ≈ domain(R)
  @test @inferred(irfft(R)) ≈ samples(xr)
  @test framerate(irfft(R)) == 1000
  @test brfft(R) ≈ 64 .* samples(xr)
  @test irfft(R, 64) ≈ samples(xr)
  @test_throws ArgumentError irfft(R, 65)
  @test irfft(fft(xr)) ≈ samples(xr)
  @test domain(R[1:2:end]) ≈ (0:16) .* (1000 / 32)
  @test irfft(R[1:2:end]) ≈ samples(xr)[1:32] .+ samples(xr)[33:64]
  @test_throws ArgumentError ifft(R)
  @test_throws ArgumentError irfft(R[1:20])
  RB = padded(R[5:10])
  RZ = zeros(ComplexF64, 33)
  RZ[5:10] .= samples(R)[5:10]
  @test irfft(RB) ≈ irfft(RZ, 64)
  @test samples(ifft(RB)) ≈ irfft(RZ, 64)
  @test RB[52] ≈ conj(RB[6])
  @test irfft(padded(fft(xr)[5:10])) ≈ irfft(RZ, 64)
  @test_throws ArgumentError padded(R[5:10], (5, 30))
  xr2 = signal(rfft(samples(xr)), 1000)
  @test irfft(xr2, 64) ≈ samples(xr)
  @test brfft(xr2, 64) ≈ 64 .* samples(xr)

  H = spectrum(samples(X), 1000)
  @test H isa SignalAnalysis.SampledSpectrum
  @test ifft(H) ≈ samples(xc)
  @test framerate(ifft(H)) == 1000
  @test irfft(spectrum(samples(R), 1000; onesided=true)) ≈ samples(xr)
  x9 = signal(randn(rng, 9), 1000)
  @test irfft(spectrum(rfft(samples(x9)), 1000; onesided=true, n=9)) ≈ samples(x9)
  @test_throws ArgumentError spectrum(samples(X), 1000; n=64)
  @test SignalAnalysis.spectruminfo(spectrum([1.0 + 0im], 1000; onesided=true)).n == 1
  @test domain(spectrum([1.0 + 0im], 1000; onesided=true)) == [0]
  @test_throws BoundsError padded(X)[2, 5]
  @test_throws ArgumentError spectrum(samples(R), 1000; onesided=true, n=66)

  xm = signal(randn(rng, ComplexF64, 16, 2), 1000)
  Xm = fft(xm)
  @test size(Xm) == (16, 2)
  @test Xm ≈ fft(samples(xm), 1)
  @test fft(xm, 1) ≈ Xm
  @test fft(xm, (1,)) ≈ Xm
  @test fft(xm, 1:1) ≈ Xm
  @test_throws ArgumentError fft(xm, 2)
  @test_throws ArgumentError fft(xm, 1:2)
  @test_throws ArgumentError fft(xm, (1, 2))
  @test ifft(Xm) ≈ samples(xm)
  @test fftshift(Xm) ≈ fftshift(samples(Xm), 1)
  @test domain(Xm) ≈ (0:15) .* (1000 / 16)
  @test Xm[:, 2] isa SignalAnalysis.SampledSpectrum
  @test domain(Xm[3:5, 1]) ≈ domain(Xm)[3:5]
  @test Xm[1:3] isa Vector{ComplexF64}
  @test Xm[2, :] isa Vector{ComplexF64}
  @test ifft(padded(Xm[3:5, :])) ≈ ifft([zeros(2, 2); samples(Xm)[3:5, :]; zeros(11, 2)], 1)
  @test padded(Xm)[0, 2] == Xm[16, 2]

  xs = signal(randn(rng, 32, 3), 1000)
  Xs = fft(xs)
  for A ∈ (xs, Xs), d ∈ (1, 2)
    for r ∈ (sum(A; dims=d), sum(abs2, A; dims=d), prod(A; dims=d), maximum(abs, A; dims=d),
             minimum(abs, A; dims=d), reduce(+, A; dims=d), mapreduce(abs, +, A; dims=d),
             count(z -> abs(z) > 1, A; dims=d), any(z -> abs(z) > 1, A; dims=d),
             first(findmax(abs, A; dims=d)), first(findmin(abs, A; dims=d)),
             mean(A; dims=d), mean(abs2, A; dims=d), var(A; dims=d), std(A; dims=d), median(abs.(A); dims=d))
      @test !(r isa SignalAnalysis.MetaArray)
    end
    @test sum(A; dims=d) ≈ sum(samples(A); dims=d)
    @test median(abs.(A); dims=d) ≈ median(abs.(samples(A)); dims=d)
  end
  @test sum(xs) ≈ sum(samples(xs))
  d = 2
  @test @inferred(sum(xs; dims=d)) isa Matrix{Float64}
  @test @inferred(mean(abs2, Xs; dims=d)) isa Matrix{Float64}

  @test_throws ArgumentError fft(xc) .* fft(signal(randn(rng, ComplexF64, 64), 2000))
  @test_throws ArgumentError X[1:4] .+ X[5:8]
  @test X[1:4] .+ samples(X[5:8]) ≈ samples(X)[1:4] .+ samples(X)[5:8]
  @test domain(X[1:4] .+ samples(X[5:8])) == domain(X[1:4])
  @test domain(samples(X[5:8]) .+ X[1:4]) == domain(X[1:4])
  @test_throws ArgumentError X[1:5] .* fft(signal(randn(rng, ComplexF64, 10), 1000))[1:5]
  @test_throws ArgumentError X[1:33] .* rfft(real(xc))
  @test_throws ArgumentError X .* xc
  @test_throws ArgumentError xc .* X
  @test_throws ArgumentError xc .+ signal(randn(rng, ComplexF64, 64), 2000)
  @test (X .* X) isa SignalAnalysis.SampledSpectrum
  @test diff(X) isa Vector{ComplexF64}
  @test diff(X) ≈ diff(samples(X))
  @test diff(Xs; dims=1) isa Matrix{ComplexF64}
  R8 = rfft(signal(randn(rng, 8), 1000))
  P8 = padded(R8)
  @test SignalAnalysis.spectruminfo(P8[2:4]).onesided
  @test SignalAnalysis.spectruminfo(P8[1:5]).onesided
  @test !SignalAnalysis.spectruminfo(P8[-2:2]).onesided
  @test !SignalAnalysis.spectruminfo(P8[5:7]).onesided
  @test domain(P8[-2:2]) ≈ (-3:1) .* 125
  @test ifft(padded(P8[-2:2])) ≈ ifft([P8[1:2]; zeros(3); P8[-2:0]])
  @test framerate(circconv(xr, xr)) == 1000

  @test fftfreq(xr) == fftfreq(64, 1000)
  @test fftfreq(256, xr) == fftfreq(256, 1000)
  @test fftfreq(xr[1:16]) == fftfreq(16, 1000)
  @test fftfreq(xm) == fftfreq(16, 1000)

  x = signal([1.0, 2.0, 3.0, 2.0, 1.0], 10.0)
  @test x.(domain(x)) == samples(x)
  @test x(domain(x)) == samples(x)
  @test x(time(1:5, x)) == samples(x)
  @test x(0.2) == x(200ms) == x(1//5) == 3.0
  @test x(-0.1) == 0 && x(0.5) == 0
  @test x(-10.0) == 0 && x(10.0) == 0
  @test x(0.45) != 0
  @test x(0.15) > 2.5
  @test x([0.15, 0.25]) == [x(0.15), x(0.25)]
  @test x(0.15, 1) == x(0.15)
  @test_throws BoundsError x(0.15, 2)
  @test_throws MethodError x(100ms:200ms)
  @test padded(x; fill=1.0)(-5.05) ≈ 1 atol=1e-3
  @test padded(x; fill=1.0)(-0.1) == 1
  y = signal([1.0 4.0; 2.0 5.0; 3.0 6.0], 10.0)
  @test y(0.1) == [2.0, 5.0]
  @test y(0.15, 2) == y(0.15)[2]
  @test y([0.0, 0.15]) == [y(0.0)'; y(0.15)']
  @test y([0.0, 0.15], 1) == y([0.0, 0.15])[:, 1]
  @test_throws BoundsError y(0.1, 3)
  x = signal(sin.(2π * 300 * (0:1999) / 1000), 1000)
  t = 0.5 .+ 0.9 .* rand(rng, 100)
  @test maximum(abs, x(t) .- sin.(2π * 300 * t)) < 2e-3
  @test x(t .* 𝓈) == x(t)
  @test x(t)[1] isa Float64
  @test signal(Float32.(samples(x)), 1000)(0.5) isa Float32
  @test signal(ComplexF32.(samples(x)), 1000)(0.5) isa ComplexF32
  @test signal(collect(1:10), 1000)(0.0015) isa Float64
  @test @inferred(x(0.5)) isa Float64
  @test @inferred(x(0.5, 1)) isa Float64
  @test @inferred(x(t)) isa Vector{Float64}
  @test @inferred(y(0.15)) isa Vector{Float64}
  @test @inferred(y(0.15, 1)) isa Float64
  @test @inferred(y([0.0, 0.15])) isa Matrix{Float64}
  @test @inferred(y([0.0, 0.15], 2)) isa Vector{Float64}

end

function test_generate()

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

function test_basic()

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

function test_dsp()

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

  gain(f, freq, fs) = abs(only(DSP.freqresp(f, [2π * inHz(freq) / fs])))
  gain(f::AbstractVector, freq, fs) = gain(PolynomialRatio(f, [1.0]), freq, fs)

  f = fir(127, 5kHz, 10kHz; fs=44.1kHz, bandstop=true)
  @test length(f) == 127
  @test gain(f, 7.5kHz, 44100) < 0.01
  @test gain(f, 0, 44100) ≈ 1 atol=0.01
  @test gain(f, 15kHz, 44100) ≈ 1 atol=0.01
  @test_throws ArgumentError fir(128, 5kHz, 10kHz; fs=44.1kHz, bandstop=true)
  @test_throws ArgumentError fir(127, 0, 5kHz; fs=44.1kHz, bandstop=true)
  @test_throws ArgumentError fir(127, 5kHz; fs=44.1kHz, bandstop=true)

  f = fir(FIRWindow(hamming(63)), 0, 5kHz; fs=44.1kHz)
  @test length(f) == 63
  @test f ≈ digitalfilter(Lowpass(10/44.1), FIRWindow(hamming(63)))
  f = fir(FIRWindow(; transitionwidth=0.05, attenuation=60), 0, 5kHz; fs=44.1kHz)
  @test gain(f, 0, 44100) ≈ 1 atol=0.01
  @test gain(f, 10kHz, 44100) < 0.01
  f = @test_logs (:warn, r"deprecated") match_mode=:any fir(127, 0, 5kHz; fs=44.1kHz, method=FIRWindow(hamming(63)))
  @test length(f) == 63

  f = iir(4, 0, 5kHz; fs=44.1kHz)
  @test f isa SecondOrderSections
  ω = range(0, π; length=64)
  @test DSP.freqresp(f, ω) ≈ DSP.freqresp(digitalfilter(Lowpass(10/44.1), Butterworth(4)), ω)
  @test gain(f, 5kHz, 44100) ≈ 1/√2 atol=1e-6
  @test gain(f, 0, 44100) ≈ 1 atol=1e-6
  @test gain(f, 15kHz, 44100) < 0.01
  f = iir(4, 15kHz; fs=44.1kHz)
  @test gain(f, 15kHz, 44100) ≈ 1/√2 atol=1e-6
  @test gain(f, 5kHz, 44100) < 0.01
  f = iir(Chebyshev1(4, 1), 5kHz, 10kHz; fs=44.1kHz)
  @test gain(f, 7.5kHz, 44100) ≈ 1 atol=0.15
  @test gain(f, 0, 44100) < 0.01
  @test gain(f, 15kHz, 44100) < 0.01
  f = iir(Elliptic(4, 0.5, 40), 5kHz, 10kHz; fs=44.1kHz, bandstop=true)
  @test gain(f, 7.5kHz, 44100) < 0.011
  @test gain(f, 0, 44100) ≈ 1 atol=0.06
  @test gain(f, 15kHz, 44100) ≈ 1 atol=0.06
  @test_throws ArgumentError iir(4, 0, 5kHz; fs=44.1kHz, bandstop=true)
  y = filt(iir(4, 0, 5kHz; fs=44.1kHz), signal(x, 44.1kHz))
  @test framerate(y) == 44100
  @test mean(y) ≈ mean(x) atol=0.1
  @test framerate(filtfilt(iir(4, 0, 5kHz; fs=44.1kHz), signal(x, 44.1kHz))) == 44100

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

  x = rand(rng, 256)
  y = rand(rng, 256)
  @test circcorr(x) ≈ circconv(conj.(circshift(reverse(x), 1)), x)
  @test circcorr(x, y) ≈ circconv(conj.(circshift(reverse(x), 1)), y)
  @test circconv(x, y) ≈ circconv(y, x)

  x = rand(rng, ComplexF64, 256)
  y = rand(rng, ComplexF64, 256)
  @test circcorr(x) ≈ circconv(conj.(circshift(reverse(x), 1)), x)
  @test circcorr(x, y) ≈ circconv(conj.(circshift(reverse(x), 1)), y)
  @test circconv(x, y) ≈ circconv(y, x)

  for j ∈ 2:10
    x = mseq(j)
    @test x isa Array{Float64}
    @test length(x) == 2^j-1
    @test all(abs.(x) .== 1)
    y = circcorr(x)
    @test y[1] ≈ length(x)
    @test all(y[2:end] .≈ -1)
    x = gmseq(j)
    @test x isa Array{ComplexF64}
    @test length(x) == 2^j-1
    @test all(abs.(x) .== 1)
    y = circcorr(x)
    @test y[1] ≈ length(x)
    @test rms(y[2:end]) ≈ 0 atol=1e-10
  end

  @test mseq((1,3)) == mseq(3)
  @test all(gmseq(3, 0) .== 1.0)

  for j ∈ [1,2,4,8]
    H = hadamard(j)
    @test H isa Matrix{Int}
    @test H' * H == 2^j * I(2^j)
    for i ∈ 0:2^j-1
      @test hadamard(i, j) == H[i+1,:]
    end
  end

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
  x3 = pll(x2, 1kHz)
  @test meanfrequency(x3[800:end]) ≈ meanfrequency(x[800:end]) atol=0.1
  e3 = √mean(abs2.(x[800:end] .- x3[800:end]))
  @test e2/e3 > 3

  t = (0:15999) / 16000
  for (a, f0) ∈ ((1.0, 1010.0), (1.0, 1100.0), (100.0, 1010.0), (0.01, 1010.0))
    x = signal(a .* cis.(2π * f0 * t .+ 0.3), 16kHz)
    y = pll(x, 1kHz)
    @test framerate(y) == 16000
    @test samples(y)[end-3999:end] ≈ cis.(2π * f0 * t[end-3999:end] .+ 0.3) atol=1e-3
  end
  x = signal(hcat(cis.(2π * 1010 * t), 5 .* cis.(2π * 990 * t)), 16kHz)
  y = pll(x, 1kHz)
  @test size(y) == size(x)
  @test samples(y)[end-3999:end, :] ≈ samples(x)[end-3999:end, :] ./ [1 5] atol=1e-3
  @test samples(pll(signal(zeros(ComplexF64, 100), 16kHz), 1kHz)) ≈ cis.(2π * 1000 * t[1:100])

  x = cw(7kHz, 1s, 44.1kHz)
  f = fir(127, 5kHz, 10kHz; fs=44.1kHz)
  x1 = filt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = filt(f, 1, x)
  @test framerate(x1) == framerate(x)
  x1 = filtfilt(f, x)
  @test framerate(x1) == framerate(x)
  x1 = filtfilt(f, [1.0], x)
  @test framerate(x1) == framerate(x)
  @test x1 ≈ filtfilt(f, x)
  x1 = filtfilt(PolynomialRatio(f, [1.0]), x)
  @test framerate(x1) == framerate(x)
  @test x1 ≈ filtfilt(f, x)
  x1 = filtfilt(f, [1.0], [x x])
  @test framerate(x1) == framerate(x)
  @test size(x1) == (nframes(x), 2)

  x = signal(cos.(2π * 10000 .* (0:19999) ./ 96000), 96000)
  x1 = filtfilt(fir(255, 5kHz, 15kHz; fs=96kHz), x)
  @test framerate(x1) == framerate(x)
  @test x1[3000:17000] ≈ x[3000:17000] atol=1e-9

  x = cw(7kHz, 1s, 44.1kHz)
  x1 = resample(x, 3//2)
  @test framerate(x1) == 3 * framerate(x) / 2
  x1 = resample(x, 3//2, [1,1,1])
  @test framerate(x1) == 3 * framerate(x) / 2
  x1 = resample([x x], 3//2)
  @test framerate(x1) == 3 * framerate(x) / 2
  x1 = resample([x x], 3//2, [1,1,1])
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

  x = randn(rng, ComplexF64, 100)
  x1 = signal(randn(rng, 1000, 2), 10kHz)
  x2 = mfilter(x, x1)
  @test isanalytic(x2)
  @test nchannels(x2) == 2

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

  z = compose(real(x), time([32.75, 64.25, 129.0], x), [1.0, -0.8, 0.6]; duration=duration(y))
  @test pow2db(energy(y - z)/ energy(y)) < -40.0
  z = compose(x, time([32.75, 64.25, 129.0], x), [1.0, -0.8, 0.6]; duration=duration(y)) |> real
  @test pow2db(energy(y - z)/ energy(y)) < -40.0
  z = compose(x, [-0.5, 0.5], [1.0, 1.0]; duration=0.3)
  @test energy(z) ≈ 0.0
  z = compose(x, [-0.05, 0.25], [1.0, 1.0]; duration=0.3)
  @test energy(z) ≈ energy(x)

  y .+= 0.1 * randn(rng, length(y))
  t, a, m = findsignal(x, y, 3; finetune=0)
  @test t ≈ [0.000775, 0.001544, 0.003125] atol=1e-5
  @test real(a) / real(a[1]) ≈ [1.0, -0.8, 0.6] atol=1e-2
  @test length(m) == 0
  t, a, m = findsignal(x, y, 3; mfo=true)
  @test t ≈ [0.000775, 0.001544, 0.003125] atol=1e-5
  @test real(a) / real(a[1]) ≈ [1.0, -0.8, 0.6] atol=1e-2
  @test length(m) > 0

  y = compose(real(x), time([32.75, 64.25, 129.0], x), [0.8, -0.7, 1.0]; duration=0.2)
  t, a, m = findsignal(x, y, 3; finetune=0)
  @test t ≈ [0.000775, 0.001544, 0.003125] atol=1e-5
  @test real(a) / real(a[3]) ≈ [0.8, -0.7, 1.0] atol=1e-2

  @test delay!([1,2,3,4], -1) == [2,3,4,0]
  @test delay!([1,2,3,4], 1) == [0,1,2,3]

  pulse = exp.(-((-8:8) ./ 2.5) .^ 2)
  x = signal(vcat(zeros(128), pulse, zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(128), [1.0, 2.0, 3.0, 2.0, 1.0], zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(128), ComplexF64.(pulse), zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(128), ComplexF64[1.0, 2.0, 3.0, 2.0, 1.0], zeros(128)), 100.0)
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(ComplexF32, 128), ComplexF32.(pulse), zeros(ComplexF32, 128)), 100.0)
  y = copy(samples(x))
  delay!(x, 1)
  @test samples(x) == circshift(y, 1)
  delay!(x, -33.0)
  @test samples(x) == circshift(y, -32)
  delay!(x, 10.7)
  delay!(x, 21.3)
  @test samples(x) ≈ y atol=1e-3

  x = signal(vcat(zeros(ComplexF32, 128), ComplexF32[1.0, 2.0, 3.0, 2.0, 1.0], zeros(ComplexF32, 128)), 100.0)
  y = copy(samples(x))
  delay!(x, 0.1s)
  delay!(x, -10)
  @test samples(x) ≈ y atol=1e-3

  x = signal(Float64[1,2,3,4,5,6,7,0],3)
  y = copy(samples(x))
  delay!(x, 0.9999999)
  delay!(x, -0.9999999)
  @test samples(x) ≈ y atol=1e-3

  x = signal(randn(rng, 500, 3), 1000)
  y = delay(x, 3.7)
  @test samples(y) ≈ x((0:499) ./ 1000 .- 0.0037)
  @test delay(x, 3.7ms) ≈ y
  @test delay(x, 3)[4:end, :] == x[1:end-3, :]
  @test all(delay(x, 3)[1:3, :] .== 0)
  @test delay(x, -3.7) ≈ x((0:499) ./ 1000 .+ 0.0037)
  @test all(delay(x, 600.5) .== 0)
  x = signal(randn(rng, Float32, 100), 10)
  @test @inferred(delay(x, 0.5)) isa typeof(x)

  x = rand(rng, ComplexF64, 4096)
  @test vec(dzt(x, length(x))) == x
  @test dzt(x, 1) ≈ fft(x) / sqrt(length(x))
  Zx = dzt(x, 64)
  @test Zx == dzt(x, 64, 64)
  @test idzt(Zx) ≈ x
  @test sum(abs2, x) ≈ sum(abs2, Zx)
  @test vec(sum(Zx; dims=1)) / sqrt(64) ≈ x[1:64]
  @test dzt(4.2x, 64) ≈ 4.2Zx
  @test abs.(dzt(circshift(x, 10), 64)) ≈ abs.(circshift(Zx, (0, 10)))
  y = rand(rng, ComplexF64, 4096)
  Zy = dzt(y, 64)
  @test vec(Zx)' * vec(Zy) ≈ x' * y
  @test dzt(4.2x + 2.7y, 64) ≈ 4.2Zx + 2.7Zy

  x = compose(mseq(12), [0.1, 0.2], [1.0, 0.7]; duration=1.0, fs=8000)
  x += 0.1 * randn(rng, size(x))
  t, a = decompose(mseq(12), x)
  @test t ≈ [0.1, 0.2]
  @test eltype(a) == Float64
  @test a ≈ [1.0, 0.7] atol=0.01
  t, a = decompose(mseq(12), x; refine=false)
  @test t ≈ [0.1, 0.2]
  @test eltype(a) == Float64
  @test a ≈ [1.0, 0.7] atol=0.01
  t, a = decompose(mseq(12), √2 * analytic(x))
  @test t ≈ [0.1, 0.2]
  @test eltype(a) == ComplexF64
  @test a ≈ [1.0, 0.7] atol=0.01
  t, a = decompose(ComplexF64.(mseq(12)), √2 * analytic(x), 2)
  @test t ≈ [0.1, 0.2]
  @test eltype(a) == ComplexF64
  @test a ≈ [1.0, 0.7] atol=0.1
  t, a, ndx = decompose(mseq(12), samples(x), 2)
  @test t ≈ [800.0, 1600.0]
  @test ndx == [801, 1601]
  @test a ≈ [1.0, 0.7] atol=0.01

end

function test_rand()

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

  for d ∈ (RedGaussian(1000), PinkGaussian(1000))
    @test length(d) == 1000
    @test length(rand(d)) == 1000
    @test size(rand(d, 3)) == (1000, 3)
    A = rand(rng, d, 3)
    @test size(A) == (1000, 3)
    @test std(A) ≈ 1 atol=2e-1
    @test A[:,1] != A[:,2]
    x = zeros(1000)
    @test Random.rand!(d, x) === x
    @test !any(iszero, x)
    @test Random.rand!(rng, d, x) === x
    @test std(x) ≈ 1 atol=2e-1
    B = zeros(1000, 3)
    @test Random.rand!(rng, d, B) === B
    @test B[:,1] != B[:,2]
  end

end

function test_array()

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

function test_tfa()

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
