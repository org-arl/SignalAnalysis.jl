# Plot recipes

Plot recipes are enabled by importing `Plots`:
```julia
using Plots
```

We provide a plot recipe to display signals:
```julia
plot(data::SampleBuf; t0=0.0, downsample=:auto, pooling=:auto, kwargs...)
```
**Optional arguments**:
- `t0=0.0`: start time (for labeling only)
- `downsample=:auto`: downsampling factor (integer or `:auto`)
- `pooling=:auto`: pooling mode (`:auto`, `:min`, `:max`, `:mean`, `:minmax`, `nothing` or function)

If the signal is too long, it is automatically downsampled in a perceptually meaningful way. The downsampling can be controlled using the `downsample` and `pooling` keywords.

We also provide convenience plot recipes for common signal processing plots:
```@docs
psd
specgram
freqresp
```
