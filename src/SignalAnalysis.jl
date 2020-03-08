module SignalAnalysis

using Requires
using DSP
using PaddedViews
using ProgressMeter

export pow2db, amp2db, db2amp, db2pow

include("basic.jl")
include("generate.jl")
include("units.jl")

function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("plot.jl")
    @require Interact="c601a237-2ae4-5e1e-952c-7a85b0c7eef1" include("iplot.jl")
  end
end

end # module
