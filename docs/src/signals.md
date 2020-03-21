# Creating & managing signals

Signals are represented by array-like data types with a series (time, frequency or space domain) per column. The `SampleBuf` data type wraps an array and carries sampling rate as metadata with it. This allows the API to be used without having to specify sampling rate at each call. However, if a user prefers to use other array-like data types, the sampling rate may be provided as a `fs` keyword argument for API calls that require sampling rate.

All code examples in this manual assume you have imported `SignalAnalysis` and `SignalAnalysis.Units`:
```julia
using SignalAnalysis, SignalAnalysis.Units
```
The `SignalAnalysis.Units` package re-exports commonly used units (`s`, `ms`, `Hz`, `kHz`, etc) from `Unitful.jl`.

### Creating and wrapping signals

Signals wrapped in a `SampleBuf` data type may be easily created using `signal()`:
```julia
x = signal(data, fs)
```
Properties such as frame rate (sampling rate), number of channels, etc may be accessed using the [`SignalBase`](https://github.com/haberdashPI/SignalBase.jl) API. The signals can be treated as arrays, but carry sampling rate metadata with them. While most operations infer metadata from the input signal, some operations may be unable to automatically infer the frame rate of the output signal. In such cases, the `@rate` and `@samerateas` macros prove handy:
```julia
# create a 1-second signal sampled at 20 kHz
y = @rate 20kHz randn(20000)

# use DSP.filt() to filter a signal and set output signal sampling rate to be
# the same as y
y = @samerateas y filt(lpf, y)    

# use DSP.resample() to resample a signal and set output signal sampling rate
# to be 2/3 that of y
y = @samerateas 2//3 y resample(y, 2//3)
```

### Processing signals block-by-block

When processing long signals, it is common to break the signal into blocks, process each block, and collect the results. To simplify this, we have `slide()` to slide a window along the data, and call a function on each window:
```@docs
slide(f::Function, s::AbstractVector, nframes, overlap=0, args...; showprogress=true)
```
The function `f` is called for each block with 3 arguments:
- `x`: windowed view of the original array
- `blknum`: block number
- `firstframe`: frame number of the first frame in the window

The results can be optionally collated by providing the return data type as the second argument of `slide()`:
```@docs
slide(f::Function, ::Type{T}, s::AbstractVector, nframes, overlap=0, args...; showprogress=true) where T
```

Blocks may optionally also be overlapped by specifying `noverlap`. Any extra arguments to `slide()` are passed on to the `f` function. An optional keyword argument `showprogress` is useful in long running operations to automatically show a progress bar.

### API reference

```@autodocs
Modules = [SignalAnalysis]
Pages   = ["signals.jl"]
```
