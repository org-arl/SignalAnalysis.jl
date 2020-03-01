module SignalAnalysis

using DSP, Plots

export pow2db, amp2db, db2amp, db2pow

include("basic.jl")
include("generate.jl")
include("plot.jl")

end # module
