module SignalAnalysis

using Requires
using DocStringExtensions

using SignalBase
using SignalBase.Units

# from SignalBase
export nframes, nchannels, sampletype, framerate, duration
export 𝓈, ms, Hz, kHz

# from DSP
export db2amp, amp2db, pow2db, db2pow, stft

# from Peaks
export findmaxima, argmaxima, peakproms, peakproms!, peakwidths, peakwidths!
export peakheights, peakheights!, filterpeaks!, findnextmaxima

const 𝓈 = Units.s

include("signals.jl")
include("basic.jl")
include("dsp.jl")
include("tfa.jl")
include("generate.jl")
include("array.jl")
include("rand.jl")

function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plot.jl")
  @require InteractiveViz="d14badfc-0adb-4d57-980e-37858d990fa5" include("iplot.jl")
end

## precompilation workload to speed up TTFX

using PrecompileTools

@setup_workload begin
  using Test
  rng = Random.MersenneTwister(0)
  include("../test/tests-core.jl")
  @compile_workload begin
    test_signals()
    test_generate()
    test_basic()
    test_dsp()
    test_rand()
    test_array()
    test_tfa()
  end
end

end # module
