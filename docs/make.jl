using Documenter
using Plots
#using InteractiveViz

push!(LOAD_PATH,"../src/")
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
      "tfa.md",
      "array.md",
      "rand.md",
      "plot.md",
      #"iplot.md"
    ]
  ]
)

deploydocs(
  repo = "github.com/org-arl/SignalAnalysis.jl.git",
  branch = "gh-pages",
  devbranch = "master",
  devurl = "dev",
  versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)
