[![CI](https://github.com/org-arl/SignalAnalysis.jl/workflows/CI/badge.svg)](https://github.com/org-arl/SignalAnalysis.jl/actions)
[![doc-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://org-arl.github.io/SignalAnalysis.jl/stable)
[![doc-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://org-arl.github.io/SignalAnalysis.jl/dev)
[![Codecov](https://codecov.io/gh/org-arl/SignalAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/org-arl/SignalAnalysis.jl)

# SignalAnalysis.jl
### Signal analysis toolbox for Julia

While a few great signal processing packages (e.g. [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl), [`SignalOperators.jl`](https://github.com/haberdashPI/SignalOperators.jl)) are available, they have limited time-frequency analysis, sonar analysis, and baseband analysis capabilities. This `SignalAnalysis.jl` package aims to fill that gap. The package has grown out of my research needs, but I hope to expand it over time to provide a wide variety of time-frequency analysis, and baseband signal processing tools.

The `SignalAnalysis.jl` works closely with, and complements other package like [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl), and [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl), to provide its functionality. While some overlap is inevitable, we mostly leverage functionality already available in those packages, and avoid duplicating it. In some cases, simpler wrappers are provided, or some symbols are re-exported for convenience. The package also works closely with [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) for common real-world units.

While this package works with most array-like data types, it uses the [`SignalBase.jl`](https://github.com/haberdashPI/SignalBase.jl) API to represent multichannel 1D signals (time, frequency or spatial domain). While the package adopts the `SampleBuf` data type from the [`SampledSignals.jl`](https://github.com/JuliaAudio/SampledSignals.jl) package to carry sampling rate information with the sampled signal, the API design allows sampling rate to be provided as a keyword argument in most cases, enabling the user to pass in any array-like data.

### Various APIs:
- [Creating & managing signals](https://org-arl.github.io/SignalAnalysis.jl/dev/signals.html)
- [Generating signals](https://org-arl.github.io/SignalAnalysis.jl/dev/generate.html)
- [Basic signal analysis](https://org-arl.github.io/SignalAnalysis.jl/dev/basic.html)
- [Signal processing](https://org-arl.github.io/SignalAnalysis.jl/dev/dsp.html)
- [Random noise generation](https://org-arl.github.io/SignalAnalysis.jl/dev/random.html)
- [Plot recipes](https://org-arl.github.io/SignalAnalysis.jl/dev/plot.html)
- [Interactive plotting](https://org-arl.github.io/SignalAnalysis.jl/dev/iplot.html)
