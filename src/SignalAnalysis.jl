module SignalAnalysis

using AxisArrays
using DSP
using RecipesBase

export pow2db, mag2db, db2amp, db2pow

include("signal.jl")
include("basic.jl")
include("generate.jl")
include("plot.jl")

end # module
