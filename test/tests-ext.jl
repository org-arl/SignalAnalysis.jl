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
  seriestypes(p) = [s[:seriestype] for s in p.series_list]
  markershapes(p) = [s[:markershape] for s in p.series_list]
  npoints(p) = [length(s[:y]) for s in p.series_list]
  y = signal(randn(rng, 100), 1000)
  p = plot(y)
  @test seriestypes(p) == [:path]
  @test markershapes(p) == [:none]
  for p in (plot(y; seriestype=:stem), plot(y; seriestype=:sticks), sticks(y))
    @test p isa Plots.Plot
    @test (:scatter, :circle) in zip(seriestypes(p), markershapes(p))
  end
  p = plot(y; seriestype=:stem, marker=:none)
  @test all(==(:none), markershapes(p))
  p = plot(y; seriestype=:stem, marker=:rect)
  @test (:scatter, :rect) in zip(seriestypes(p), markershapes(p))
  p = scatter(y)
  @test seriestypes(p) == [:scatter]
  @test npoints(p) == [100]
  p = scatter([y -y])
  @test seriestypes(p) == [:scatter, :scatter]
  y = signal(randn(rng, 20000), 1000)
  @test npoints(plot(y)) == [10000]
  @test npoints(scatter(y)) == [5000]
  @test npoints(scatter(y; pooling=:minmax)) == [10000]
  @test npoints(scatter(y; pooling=:mean)) == [5000]
  @test npoints(scatter(y; downsample=1)) == [20000]
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
  p = freqz(f)
  @test p isa Plots.Plot
  @test length(p.subplots) == 2
  p = freqz(f; xscale=:log10)
  @test p isa Plots.Plot
  p = freqz(f, [1])
  @test p isa Plots.Plot
  p = freqz(f; fs=44.1kHz)
  @test p isa Plots.Plot
  maxdb = maximum(p[1][1][:y])
  @test all(Plots.ylims(p[1]) .≈ (maxdb - 30, maxdb + 3))
  p = freqz(f; dbrange=nothing)
  @test p isa Plots.Plot
  p = freqz(f; show=:mag)
  @test length(p.subplots) == 1
  @test p[1][:yaxis][:guide] == "Response (dB)"
  p = freqz(f; show=:phase)
  @test length(p.subplots) == 1
  @test p[1][:yaxis][:guide] == "Phase (radians)"
  @test all(abs.(diff(p[1][1][:y])) .< π)
  p = freqz(f; show=:phase, unwrap=false)
  @test all(-π .<= p[1][1][:y] .<= π)
  @test_throws ArgumentError freqz(f; show=:foo)
  p = freqz(iir(4, 0, 5kHz; fs=44.1kHz); fs=44.1kHz, ylims=[(-60, 5) (-10, 0)])
  @test Plots.ylims(p[1]) == (-60, 5)
  @test Plots.ylims(p[2]) == (-10, 0)
  p = freqz(iir(Elliptic(4, 0.5, 40), 5kHz, 10kHz; fs=44.1kHz, bandstop=true); fs=44.1kHz)
  @test p isa Plots.Plot
  p = @test_logs (:warn, r"deprecated") plotfreqresp(f)
  @test p isa Plots.Plot

end
