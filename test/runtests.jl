using Test, Statistics, LinearAlgebra, DSP, DSP.Util, FFTW
using Plots
using WAV
using SignalAnalysis
using SignalAnalysis.Units
using StableRNGs

# used by core tests for repeatable results
rng = StableRNG(0)

# core tests
include("tests-core.jl")

@testset "signals" begin
  test_signals()
end

@testset "generate" begin
  test_generate()
end

@testset "basic" begin
  test_basic()
end

@testset "dsp" begin
  test_dsp()
end

@testset "rand" begin
  test_rand()
end

@testset "array" begin
  test_array()
end

@testset "tfa" begin
  test_tfa()
end

# WAV extension tests
@testset "wav" begin
  mktempdir() do dir
    filename = joinpath(dir, "test.wav")
    x = signal(randn(rng, 8000), 8000)
    wavwrite(samples(x), filename; Fs=8000, nbits=64)
    y = signal(filename)
    @test framerate(y) == 8000
    @test nframes(y) == 8000
    @test vec(samples(y)) ≈ samples(x)
    y = signal(filename; start=101, nsamples=1000)
    @test nframes(y) == 1000
    @test vec(samples(y)) ≈ samples(x)[101:1100]
    y = signal(filename; start=101)
    @test nframes(y) == 7900
  end
end

# plotting extension tests
include("tests-ext.jl")
