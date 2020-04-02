module SignalAnalysis

using Requires
using DocStringExtensions

using SignalBase
using SignalBase.Units
export nframes, nchannels, sampletype, framerate, duration

include("signals.jl")

include("basic.jl")
include("dsp.jl")
include("generate.jl")
include("rand.jl")

function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("plot.jl")
    include("iplot.jl")
  end
end

end # module
