module SignalAnalysis

using Requires
using Reexport
using DocStringExtensions

@reexport using SignalBase
@reexport using SignalBase.Units

@reexport using DSP
@reexport using FFTW
@reexport using Peaks
@reexport using Statistics
@reexport using LinearAlgebra

export 𝓈, ms, Hz, kHz, °

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

@compile_workload begin
  include("precompile.jl")
end

end # module
