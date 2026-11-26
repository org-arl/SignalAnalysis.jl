module SignalAnalysis

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
include("spectrum.jl")
include("basic.jl")
include("dsp.jl")
include("tfa.jl")
include("generate.jl")
include("array.jl")
include("rand.jl")

# implemented in the Plots extension (see ext/PlotsExt.jl)
export psd, psd!, specgram, specgram!, freqz, freqz!, plotfreqresp, plotfreqresp!

"""
Plots the power spectral density of data. Requires `Plots` to be loaded.
"""
function psd end
function psd! end

"""
Plots a spectrogram of the data. Requires `Plots` to be loaded.
"""
function specgram end
function specgram! end

"""
Plots frequency response of a digital filter. Requires `Plots` to be loaded.
"""
function freqz end
function freqz! end

function plotfreqresp(args...; kwargs...)
  @warn "`plotfreqresp` is deprecated, use `freqz` instead" maxlog=1
  freqz(args...; kwargs...)
end

function plotfreqresp!(args...; kwargs...)
  @warn "`plotfreqresp!` is deprecated, use `freqz!` instead" maxlog=1
  freqz!(args...; kwargs...)
end

# implemented in the InteractiveViz extension (see ext/InteractiveVizExt.jl)
export iplot, iplot!, ispecgram

"""
Plots interactive timeseries of the signal. Requires `InteractiveViz` to be loaded.
"""
function iplot end

"""
Plots interactive timeseries of the signal over a previous plot. Requires
`InteractiveViz` to be loaded.
"""
function iplot! end

"""
Plots interactive spectrogram of the signal. Requires `InteractiveViz` to be loaded.
"""
function ispecgram end

## precompilation workload to speed up TTFX

using PrecompileTools

@compile_workload begin
  include("precompile.jl")
end

end # module
