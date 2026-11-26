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
A number without units takes the units of the other end of the range, so `x[250:500ms]` is the segment from 250 to 500 milliseconds. A time range may also be given as a tuple, `x[(250ms, 500ms)]`, in which case both ends must have units, since a tuple of plain numbers would look like indices.

### Creating and wrapping signals

Signals wrapped in a `SampledSignal` data type may be easily created using `signal()`:
```julia
x = signal(data, fs)
```
Properties such as frame rate (sampling rate), number of channels, etc may be accessed using the [`SignalBase`](https://github.com/haberdashPI/SignalBase.jl) API. The signals can be treated as arrays, but carry sampling rate metadata with them. Most operations infer the sampling rate of their output from the input signal, and combining signals with different sampling rates throws an error. Reductions with `dims` (e.g. `sum(x; dims=2)` to combine channels) return plain arrays; use `samerateas(x, sum(x; dims=2))` to keep the sampling rate. Common `DSP.jl` functions retain the sampling rate when applied to signals, or adjust it when resampling:
```julia
# create a 1-second signal sampled at 20 kHz
y = signal(randn(20000), 20kHz)

# design a lowpass filter, and filter the signal, retaining its sampling rate
lpf = fir(127, 0, 5kHz; fs=20kHz)
y1 = filt(lpf, y)

# zero-phase filtering also retains the sampling rate
y1 = filtfilt(lpf, y)

# resampling adjusts the sampling rate, to 13333.33 Hz here
y1 = resample(y, 2//3)
```
The FFT of a signal is a spectrum, which also carries sampling rate information, as described in [Spectra](@ref).

### Padding signals

`padded()` pads a signal with zeros (or another `fill` value). Without an amount of padding, the signal is padded indefinitely in both directions. The result is a view with the same length and indices as the signal, but any index is valid, and indices outside the signal give the fill value. This is the simplest way to access a signal beyond its ends, e.g. in loops that access `y[i-k]`:
```julia
x = signal(randn(20000), 20kHz)
y = padded(x)               # zeros before and after x, without limit
y[-1000]                    # 0.0
y[1]                        # same as x[1]
y = padded(x; delay=10)     # x delayed by 10 samples, so y[i] == x[i-10]
```
With a given amount of padding, a signal is extended by that many samples. The result is a new signal, indexed from 1, with the padding prepended and appended. Padding on the left therefore delays the signal:
```julia
y = padded(x, (100, 200))   # 100 zeros before x and 200 zeros after it
y[101]                      # same as x[1], now at time 100/20000 = 5 ms
```
For plain arrays (rather than signals), finite padding instead gives a view that keeps the original indices, so `padded(v, 100)[-99]` is the first sample of the padding and `padded(v, 100)[1] == v[1]`.

### Values between samples

A signal may be called as a function of time, to treat it as a continuous-time signal. The time is in seconds, or given with units, and time 0 is the first sample:
```julia
x = signal(randn(20000), 20kHz)
x(0.1)                      # value at 100 ms, same as x[2001]
x(100.03ms)                 # value between two samples
x(0:1e-5:0.01)              # values at many times, as a vector
```
Values between samples are computed by band-limited interpolation, using a Kaiser-windowed sinc kernel spanning about 37 samples. This is the same design as the default filter of `resample()`. Frequencies above about 90% of the Nyquist frequency are attenuated, so signals that are to be evaluated between samples should have little content in that band. At sample times, the samples are returned unchanged. Samples beyond the ends of the signal are taken to be zero (or the fill value of a [`padded`](@ref) signal), so values within about 18 samples outside the signal may be non-zero.

For a multichannel signal, `x(t)` returns a vector with one value per channel, and `x(t, ch)` returns the value of channel `ch`. A vector of times gives a matrix with one row per time and one column per channel.

Each value takes about 37 multiply-adds per channel to compute. To evaluate a whole signal on a new sampling grid, `resample()` is much faster. To delay a signal by a fractional number of samples, use [`delay`](@ref), which uses the same interpolation. A signal delayed by a time `τ` (e.g. `delay(x, 0.25ms)`) has the value `x(t - τ)` at each sample time `t`. Note that a delay given as a plain number is in samples (e.g. `delay(x, 2.5)`), while a time given to `x(t)` as a plain number is in seconds.

### API reference

```@autodocs
Modules = [SignalAnalysis]
Pages   = ["signals.jl"]
```
