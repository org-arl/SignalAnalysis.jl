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
