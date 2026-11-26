# Plot recipes

Plot recipes are enabled by importing `Plots`:
```julia
using Plots
```

We provide a plot recipe to display signals:
```julia
plot(data::SampledSignal; t0=0.0, downsample=:auto, pooling=:auto, kwargs...)
```
**Optional arguments**:
- `t0=0.0`: start time (for labeling only)
- `downsample=:auto`: downsampling factor (integer or `:auto`)
- `pooling=:auto`: pooling mode (`:auto`, `:min`, `:max`, `:mean`, `:minmax`, `nothing` or function)
- `seriestype=:line`: series type, e.g. `:line`, `:stem` or `:scatter`

If the signal is too long, it is automatically downsampled in a perceptually meaningful way. The downsampling can be controlled using the `downsample` and `pooling` keywords. By default, line plots use `:minmax` pooling, and other series types use no pooling.

A signal is plotted as a line by default. Stem and scatter plots show the individual samples:
```julia
plot(x; seriestype=:stem)
scatter(x)
```
Stem plots have circle markers by default, and `marker=:none` removes them.

Spectra (see [Spectra](@ref)) are plotted as magnitude (in dB) against frequency:
```julia
plot(X::SampledSpectrum; db=true, wrap=true, kwargs...)
```
**Optional arguments**:
- `db=true`: plot magnitude in dB, or linear magnitude if `false`
- `wrap=true`: show frequencies wrapped to `[-fs/2, fs/2)`, or unwrapped (`0` to `fs`) if `false`
- `ylims=50`: a number gives the dB range below the maximum; limits may also be given as usual

For spectra, power spectral densities and filter frequency responses, a single number for `ylims` gives the dB range of the plot below its maximum, while a tuple gives the limits as usual.

We also provide convenience plot recipes for common signal processing plots:
```@docs
psd
specgram
freqz
```
