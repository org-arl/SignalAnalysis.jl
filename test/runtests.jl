using Test, Statistics, LinearAlgebra, DSP, DSP.Util, FFTW
using Plots
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

# plotting extension tests
include("tests-ext.jl")
