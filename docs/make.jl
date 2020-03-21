using Documenter, Plots

push!(LOAD_PATH,"../../src/")
using SignalAnalysis

makedocs(
  sitename = "SignalAnalysis.jl",
  format = Documenter.HTML(prettyurls = false),
  linkcheck = !("skiplinks" in ARGS),
  pages = Any[
    "Home" => "index.md",
    "Manual" => Any[
      "signals.md",
      "generate.md",
      "basic.md",
      "dsp.md",
      "rand.md",
      "plot.md",
      "iplot.md"
    ]
  ]
)

deploydocs(
  repo = "github.com/org-arl/SignalAnalysis.jl.git",
)
