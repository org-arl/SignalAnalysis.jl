# Spectra

The FFT of a signal is a _spectrum_: an array of FFT coefficients that also knows the sampling rate of the signal it came from, the length of the FFT, and which bins of the FFT grid it holds. With this information, a spectrum knows the frequency of each of its bins, and can be transformed back into a signal at the right sampling rate, even after it has been sliced, shifted, concatenated or padded.

### From signals to spectra

`fft` and `rfft` of a signal return spectra:
```julia
x = signal(randn(8000), 8kHz)   # 1 second of noise sampled at 8 kHz
X = fft(x)                      # two-sided spectrum with 8000 bins
R = rfft(x)                     # one-sided spectrum with 4001 bins
framerate(X)                    # 8000.0, the sampling rate of x
samples(X)                      # the underlying FFT coefficients
```
Each channel of a multichannel signal is transformed separately, along time. `fft(x, 1)` is accepted as a synonym for `fft(x)`, but other regions are not, since they would transform across channels. To transform along other dimensions, use the underlying samples, e.g. `fft(samples(x), 2)`.

A spectrum is an array, so it can be used in computations like any other array. Combining a spectrum with plain arrays or numbers, or with another spectrum covering the same bins of the same FFT grid, gives a spectrum. Combining spectra with different sampling rates, FFT lengths or bins, or a spectrum with a signal, throws an error, since this is usually a mistake. To combine different bands deliberately, use `samples()` on one of them, so that the result keeps the bins of the other:
```julia
X .* conj.(Y)                       # spectrum, if X and Y cover the same bins
X[1kHz:1999Hz] .+ samples(X[2kHz:2999Hz])   # spectrum with bins from 1 kHz to 1999 Hz
```
Reductions with `dims` (e.g. `sum(abs2, X; dims=1)` or `mean(X; dims=2)`) return plain arrays, as do functions such as `diff` whose results no longer correspond to the bins of the spectrum.

### Frequencies

`domain` returns the frequency of each bin, as a range. The frequencies are not wrapped, so an unshifted two-sided spectrum spans `0` to `fs`, with the negative frequencies in the upper half. `fftfreq` returns the same frequencies wrapped to the interval `[-fs/2, fs/2)`:
```julia
domain(X)                # 0.0:1.0:7999.0 Hz
fftfreq(X)               # 0, 1, …, 3999, -4000, …, -1 Hz
```
`freq` converts a bin index to its frequency, and `tobin` converts a frequency to the index of the nearest bin, like `time` and `toframe` do for signals:
```julia
freq(1001, X)            # 1000.0 Hz
tobin(1kHz, X)           # 1001
tobin((1kHz, 2kHz), X)   # 1001:2001
tobin(1kHz, 8kHz, 8000)  # 1001, for an 8000-point FFT at 8 kHz
```

### Values between bins

A spectrum may be called as a function of frequency, to get its value between bins. The frequency is in Hz, or given with units:
```julia
X(1000.5)                # value halfway between two bins
X(1.0005kHz)             # same value
X(990:0.1:1010)          # values at many frequencies, as a vector
```
The value is the discrete-time Fourier transform of the signal that the spectrum represents, so it is exact, and matches the bins at bin frequencies. This is useful, for example, to find the frequency of a spectral peak more precisely than the bin spacing. The spectrum must be complete: all its bins (shifted or not), all non-negative frequency bins of a one-sided spectrum, or a padded spectrum. Each value takes time proportional to the FFT length, so for values on a fine regular grid, zero-padding the signal before the FFT is faster. For a multichannel spectrum, `X(f)` returns one value per channel, and `X(f, ch)` the value of channel `ch`.

### Selecting bands

Spectra can be indexed by frequency, in the same way that signals can be indexed by time. The result is a spectrum that still knows the frequencies of its bins:
```julia
B = X[1kHz:2kHz]         # bins from 1 kHz to 2 kHz
B = X[1:2kHz]            # the same, since a number without units takes the units of the other end
domain(B)                # 1000.0:1.0:2000.0 Hz
X[1kHz:2kHz] = zeros(1001)   # zero the band from 1 kHz to 2 kHz
```
In an unshifted two-sided spectrum, negative frequencies are stored in the upper half of the spectrum, so a band that includes negative frequencies is not a contiguous block of bins, and cannot be selected. Shift the spectrum first, or use a padded spectrum (see below):
```julia
fftshift(X)[-1kHz:1kHz]  # bins from -1 kHz to 1 kHz
padded(X)[-1kHz:1kHz]    # the same, without reordering X
```
`fftshift` and `ifftshift` move a complete two-sided spectrum between the unshifted grid (`0` to `fs`) and the centred one (`-fs/2` to `fs/2`).

Slicing a spectrum with a range of indices keeps track of its frequencies. A stride `s` is allowed if it divides both the FFT length and the index of the first bin. The result is the spectrum of an FFT that is `s` times shorter, i.e. of the signal folded (time-aliased) to `1/s` of its length:
```julia
X[1:2:end]               # bins 0, 2, 4, … Hz: the spectrum of a 4000-point FFT
ifft(X[1:2:end])         # x[1:4000] + x[4001:8000]
```
Other strides throw an error. Selecting bins with an index vector or a logical mask returns a plain array, since the bins no longer form a regular grid. Selecting channels (e.g. `X[:, 2]`) keeps the spectrum.

### Back to time

`ifft` and `irfft` transform a spectrum back into a signal at the original sampling rate. `irfft` uses the bins with non-negative frequencies, so it accepts both one-sided and two-sided spectra:
```julia
y = ifft(X)              # complex signal at 8 kHz
y = irfft(R)             # real signal at 8 kHz
y = irfft(X)             # the same, from the two-sided spectrum
```
An inverse transform needs a complete spectrum: all the bins of an unshifted FFT grid. A slice or a shifted spectrum must be padded first. A padded spectrum is a view of the spectrum, with unchanged indexing, that extends it indefinitely: with zeros within one period of the FFT grid, and repeating periodically beyond it, as the frequencies of an FFT do. A one-sided spectrum is extended by conjugate symmetry, as for the spectrum of a real signal. An inverse transform of a padded spectrum selects one period of it:
```julia
y = irfft(padded(X[1kHz:2kHz]))   # real signal band-limited to 1–2 kHz
y = irfft(padded(R[1kHz:2kHz]))   # the same, from the one-sided spectrum
y = ifft(padded(X[1kHz:2kHz]))    # complex signal with only positive frequencies from 1–2 kHz
y = ifft(padded(fftshift(X)))     # no need for ifftshift
```
A spectrum may also be padded with a given number of zero bins on each side, e.g. `padded(B, (10, 20))`, as long as it stays within one period of the FFT grid. The result is a new spectrum indexed from 1, whose frequencies are adjusted for the added bins.

### Building spectra

`spectrum` creates a spectrum from FFT coefficients, for example to compute the impulse response of a filter specified in the frequency domain:
```julia
f = fftfreq(256, 8000)                          # frequencies of a 256-point FFT at 8 kHz
H = spectrum(Float64.(abs.(f) .≤ 1000), 8kHz)   # ideal lowpass filter with 1 kHz cutoff
h = ifft(H)                                     # its (circular) impulse response at 8 kHz
```
One-sided spectra are created with `spectrum(H, fs; onesided=true)`, adding `n` for an odd FFT length.

Spectra may be concatenated. If all the pieces are spectra from the same FFT grid, and they form contiguous segments of it (for `vcat`) or have identical bins (for `hcat`), the result is a spectrum. Pieces that do not fit together throw an error, and if any piece is a plain array, the result is a plain array:
```julia
Y = [X[1:1000]; X[1001:8000]]   # a spectrum again, equal to X
```

### Plotting

When `Plots` is loaded, spectra can be plotted directly, showing magnitude (in dB) against frequency:
```julia
plot(X)                  # frequencies wrapped to -fs/2 … fs/2
plot(X; wrap=false)      # unwrapped frequencies, 0 … fs
plot(X; db=false)        # linear magnitude
plot(X; ylims=30)        # 30 dB below the maximum (default 50 dB)
```

### API reference

```@autodocs
Modules = [SignalAnalysis]
Pages   = ["spectrum.jl"]
```
