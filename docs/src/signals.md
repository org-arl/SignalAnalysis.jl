# Creating & managing signals

Signals are represented by array-like data types with a series (time, frequency or space domain) per column. The `SampledSignal` data type wraps an array and carries sampling rate as metadata with it. This allows the API to be used without having to specify sampling rate at each call. However, if a user prefers to use other array-like data types, the sampling rate may be provided as a `fs` keyword argument for API calls that require sampling rate.

All code examples in this manual assume you have imported `SignalAnalysis` and `SignalAnalysis.Units`:
```julia
using SignalAnalysis, SignalAnalysis.Units
```
The `SignalAnalysis.Units` package re-exports commonly used units (`s`, `ms`, `Hz`, `kHz`, etc) from `Unitful.jl`. In addition, a variable `𝓈` is always exported and is an alias for the `s` unit from `SignalAnalysis.Units`. This allows indexing of signals using time in seconds:
```julia
# create a 1-second signal sampled at 20 kHz
x = signal(randn(20000), 20kHz)

# get a signal segment from 0.25 to 0.5 seconds:
y = x[0.25:0.5𝓈]
```

### Creating and wrapping signals

Signals wrapped in a `SampledSignal` data type may be easily created using `signal()`:
```julia
x = signal(data, fs)
```
Properties such as frame rate (sampling rate), number of channels, etc may be accessed using the [`SignalBase`](https://github.com/haberdashPI/SignalBase.jl) API. The signals can be treated as arrays, but carry sampling rate metadata with them. While most operations infer metadata from the input signal, some operations may be unable to automatically infer the frame rate of the output signal. We provide some handy wrappers around common `DSP.jl` functions to aid with rate inference:
```julia
# create a 1-second signal sampled at 20 kHz
y = signal(randn(20000), 20kHz)

# use DSP.filt() to filter a signal but retainin sampling rate
y = sfilt(lpf, y)

# use DSP.filtfilt() to filter a signal but retainin sampling rate
y = sfiltfilt(lpf, y)

# use DSP.resample() to resample a signal and infer sampling rate
y = sresample(y, 2//3)
```

### API reference

```@autodocs
Modules = [SignalAnalysis]
Pages   = ["signals.jl"]
```
