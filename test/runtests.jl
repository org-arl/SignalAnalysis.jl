using Test, Statistics, LinearAlgebra, DSP, DSP.Util
using Plots
using SignalAnalysis
using SignalAnalysis.Units
using StableRNGs

# used by core tests for repeatable results
rng = StableRNG(0)

# core tests
include("tests-core.jl")
@testset "signals" test_signals()
@testset "generate" test_generate()
@testset "basic" test_basic()
@testset "dsp" test_dsp()
@testset "rand" test_rand()
@testset "array" test_array()
@testset "tfa" test_tfa()

# plotting extension tests
include("tests-ext.jl")
