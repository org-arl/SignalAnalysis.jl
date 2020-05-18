# SignalAnalysis.jl
### Signal analysis toolbox for Julia

```@meta
CurrentModule = SignalAnalysis
```

While a few great signal processing packages (e.g. [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl), [`SignalOperators.jl`](https://github.com/haberdashPI/SignalOperators.jl)) are available, they have limited time-frequency analysis, sonar analysis, and baseband analysis capabilities. This `SignalAnalysis.jl` package aims to fill that gap. The package has grown out of my research needs, but I hope to expand it over time to provide a wide variety of time-frequency analysis, and baseband signal processing tools.

While this package works with most array-like data types, it uses the [`SignalBase.jl`](https://github.com/haberdashPI/SignalBase.jl) API to represent multichannel 1D signals (time, frequency or spatial domain). While the package adopts the `SampleBuf` data type from the [`SampledSignals.jl`](https://github.com/JuliaAudio/SampledSignals.jl) package to carry sampling rate information with the sampled signal, the API design allows sampling rate to be provided as a keyword argument in most cases, enabling the user to pass in any array-like data.

### APIs

- [Creating & managing signals](@ref)
- [Generating signals](@ref)
- [Basic signal analysis](@ref)
- [Signal processing](@ref)
- [Array signal processing](@ref)
- [Random noise generation](@ref)
- [Plot recipes](@ref)
- [Interactive plotting](@ref)

### Quick links

- [`SignalAnalysis.jl`](https://github.com/org-arl/SignalAnalysis.jl)
- [`SignalBase.jl`](https://github.com/haberdashPI/SignalBase.jl)
- [`SampledSignals.jl`](https://github.com/JuliaAudio/SampledSignals.jl)
- [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl)
- [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl)
- [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
