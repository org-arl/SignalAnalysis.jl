export fir, iir, removedc, removedc!, demon
export upconvert, downconvert, rrcosfir, rcosfir
export mseq, gmseq, circconv, circcorr, goertzel, pll, hadamard
export mfilter, findsignal, dzt, idzt
export istft, whiten, filt, filtfilt, resample, delay, delay!, compose, decompose

function responsetype(f1, f2, fs, bandstop)
  fs = inHz(fs)
  if f1 == 0
    bandstop && throw(ArgumentError("bandstop filters need both band edges f1 and f2"))
    Lowpass(2 * inHz(f2) / fs)
  elseif f2 === nothing || inHz(f2) == fs/2
    bandstop && throw(ArgumentError("bandstop filters need both band edges f1 and f2"))
    Highpass(2 * inHz(f1) / fs)
  elseif bandstop
    Bandstop(2 * inHz(f1) / fs, 2 * inHz(f2) / fs)
  else
    Bandpass(2 * inHz(f1) / fs, 2 * inHz(f2) / fs)
  end
end

"""
    fir(n, f1, f2=nothing; fs=2.0, bandstop=false)
    fir(method, f1, f2=nothing; fs=2.0, bandstop=false)

Designs a FIR filter with a passband from `f1` to `f2`. If frame rate `fs` is not
specified, `f1` and `f2` are given in normalized units (1.0 being Nyquist). If `f1`
is 0, the designed filter is a lowpass filter, and if `f2` is `nothing` then it is
a highpass filter. If `bandstop` is `true`, the band from `f1` to `f2` is rejected
instead.

The filter is returned as a vector of taps, and can be applied using `filt(f, x)`
or `filtfilt(f, x)`. If an integer `n` is given, a `n`-tap Hann-windowed filter is
designed. Otherwise
`method` is a `DSP.FIRWindow` design method, such as `FIRWindow(hamming(63))` or
`FIRWindow(; transitionwidth=0.05, attenuation=60)`. Highpass and bandstop filters
need an odd number of taps.

This method is a convenience wrapper around `DSP.digitalfilter`.

# Examples:
```julia-repl
julia> lpf = fir(127, 0, 10kHz; fs=44.1kHz)   # design a lowpass filter
127-element Array{Float64,1}:
  ⋮

julia> hpf = fir(127, 10kHz; fs=44.1kHz)      # design a highpass filter
127-element Array{Float64,1}:
  ⋮

julia> bpf = fir(127, 1kHz, 5kHz; fs=44.1kHz) # design a bandpass filter
127-element Array{Float64,1}:
  ⋮

julia> bsf = fir(127, 1kHz, 5kHz; fs=44.1kHz, bandstop=true) # design a bandstop filter
127-element Array{Float64,1}:
  ⋮

julia> lpf = fir(FIRWindow(hamming(63)), 0, 10kHz; fs=44.1kHz) # use a Hamming window
63-element Array{Float64,1}:
  ⋮
```
"""
function fir(n::Integer, f1, f2=nothing; fs=2.0, bandstop=false, method=nothing)
  if method !== nothing
    @warn "`fir(n, f1, f2; method)` is deprecated, use `fir(method, f1, f2)` instead" maxlog=1
    return fir(method, f1, f2; fs, bandstop)
  end
  fir(FIRWindow(hanning(n)), f1, f2; fs, bandstop)
end

fir(method::FIRWindow, f1, f2=nothing; fs=2.0, bandstop=false) =
  digitalfilter(responsetype(f1, f2, fs, bandstop), method)

"""
    iir(n, f1, f2=nothing; fs=2.0, bandstop=false)
    iir(method, f1, f2=nothing; fs=2.0, bandstop=false)

Designs an IIR filter with a passband from `f1` to `f2`. If frame rate `fs` is not
specified, `f1` and `f2` are given in normalized units (1.0 being Nyquist). If `f1`
is 0, the designed filter is a lowpass filter, and if `f2` is `nothing` then it is
a highpass filter. If `bandstop` is `true`, the band from `f1` to `f2` is rejected
instead.

If an integer `n` is given, an order `n` Butterworth filter is designed. Otherwise
`method` is a DSP.jl IIR design method, such as `Butterworth(n)`,
`Chebyshev1(n, ripple)`, `Chebyshev2(n, ripple)` or `Elliptic(n, rp, rs)`.
Bandpass and bandstop filters have twice the order of the design method.

The filter is returned as `DSP.SecondOrderSections`, for numerical stability, and
can be applied using `filt(f, x)` or `filtfilt(f, x)`. Numerator and denominator
coefficients can be obtained using `convert(PolynomialRatio, f)`, but filtering with
them may be inaccurate or unstable for high order or narrowband filters.
This method is a convenience wrapper around `DSP.digitalfilter`.

# Examples:
```julia-repl
julia> lpf = iir(4, 0, 1kHz; fs=8kHz)                       # Butterworth lowpass filter
julia> hpf = iir(4, 100Hz; fs=8kHz)                         # Butterworth highpass filter
julia> bpf = iir(Chebyshev1(4, 1), 1kHz, 2kHz; fs=8kHz)     # Chebyshev bandpass filter
julia> bsf = iir(Elliptic(4, 0.5, 40), 950Hz, 1050Hz; fs=8kHz, bandstop=true)
julia> y = filt(lpf, x)                                     # apply the filter
julia> freqz(bpf; fs=8kHz)                                  # plot the frequency response
```
"""
iir(n::Integer, f1, f2=nothing; fs=2.0, bandstop=false) =
  iir(Butterworth(n), f1, f2; fs, bandstop)

iir(method::FilterCoefficients{:s}, f1, f2=nothing; fs=2.0, bandstop=false) =
  convert(SecondOrderSections, digitalfilter(responsetype(f1, f2, fs, bandstop), method))

"""
$(SIGNATURES)
DC removal filter. Parameter `α` controls the cutoff frequency. Implementation
based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).

See also: [`removedc`](@ref)
"""
function removedc!(s; α=0.95)
  for k = 1:size(s,2)
    for j = 2:size(s,1)
      s[j,k] += α*s[j-1,k]
    end
    s[2:end,k] .+= -s[1:end-1,k]
  end
  s *= sqrt(α)
  return s
end

"""
$(SIGNATURES)
DC removal filter. Parameter `α` controls the cutoff frequency. Implementation
based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).

See also: [`removedc!`](@ref)
"""
removedc(s; α=0.95) = removedc!(copy(s); α=α)

"""
$(SIGNATURES)
Estimates DEMON spectrum. The output is highpass filtered with a `cutoff`
frequency and downsampled. Supported downsampling methods are `:rms` (default),
`:mean` and `:fir`.
"""
function demon(x; fs=framerate(x), downsample=250, method=:rms, cutoff=1.0)
  local y
  fs /= downsample
  for k = 1:size(x,2)
    if downsample == 1
      y1 = abs.(hilbert(x[:,k]))
    elseif method == :rms
      y1 = sqrt.(mean.(Periodograms.arraysplit(abs2.(hilbert(x[:,k])), downsample, 0)))
    elseif method == :mean
      y1 = mean.(Periodograms.arraysplit(abs.(hilbert(x[:,k])), downsample, 0))
    elseif method == :fir
      aaf = fir(127, 0, 0.48fs; fs=fs)
      y1 = filtfilt(aaf, abs.(hilbert(x[:,k])))[1:downsample:end]
    else
      throw(ArgumentError("Unknown method"))
    end
    if k == 1
      y = zeros(length(y1), size(x,2))
    end
    y[:,k] .= y1
  end
  maxflen = length(y)÷4
  mod(maxflen, 2) == 0 && (maxflen += 1)
  hpf = fir(min(127, maxflen), cutoff; fs=fs)
  signal(filtfilt(hpf, y), fs)
end

function demon(x::AbstractVector{T}; fs=framerate(x), downsample=250, method=:rms, cutoff=1.0) where T
  y = @view samples(x)[:,1:1]
  z = demon(y, fs=fs, downsample=downsample, method=method, cutoff=cutoff)
  signal(dropdims(samples(z), dims=2), framerate(z))
end

"""
$(SIGNATURES)
Converts baseband signal with `sps` symbols per passband sample to a real
passband signal centered around carrier frequency `fc`.
"""
function upconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  if sps != 1
    pad = cld(length(pulseshape), 2*sps) - 1
    s = vcat(zeros(pad), complex.(s), zeros(pad))
    s = signal(resample(s, sps, pulseshape), sps*fs)
  end
  √2 * real.(s .* cis.(2π * inHz(fc) * domain(s)))
end

function upconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  mapreduce(hcat, eachcol(s)) do x
    upconvert(x, sps, fc, pulseshape; fs=fs)
  end[:,:]
end

"""
$(SIGNATURES)
Converts passband signal centered around carrier frequency `fc` to baseband,
and downsamples it by a factor of `sps`. If the `pulseshape` is specified to
be `nothing`, downsampling is performed without filtering.
"""
function downconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  s = signal(analytic(s), fs)
  s = s .* cis.(-2π * inHz(fc) * domain(s))
  sps == 1 && return signal(s, fs)
  pulseshape == nothing && return signal(s[1:sps:end,:], fs/sps)
  signal(resample(s, 1//sps, pulseshape), fs/sps)
end

function downconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  mapreduce(hcat, eachcol(s)) do x
    downconvert(x, sps, fc, pulseshape; fs=fs)
  end[:,:]
end

"""
$(SIGNATURES)
Root-raised cosine filter.
"""
function rrcosfir(β, sps, span = β < 0.68 ? 33-floor(Int, 44β) : 4)
  # default span based on http://www.commsys.isy.liu.se/TSKS04/lectures/3/MichaelZoltowski_SquareRootRaisedCosine.pdf
  delay = fld(span*sps, 2)
  t = collect(-delay:delay)/sps
  h = Array{Float64}(undef, size(t))
  for i ∈ 1:length(t)
    if t[i] == 0
      h[i] = (1 + β*(4/π - 1))/sps
    elseif abs(t[i]) == 1/(4β)
      h[i] = β/(√2*sps) * ((1+2/pi)*sin(π/(4β)) + (1-2/pi)*cos(π/(4β)))
    else
      h[i] = (sin(π*t[i]*(1-β)) + 4β*t[i]*cos(π*t[i]*(1+β))) / (π*t[i]*(1 - (4β*t[i])^2)) / sps
    end
  end
  h / √sum(h.^2)
end

"""
$(SIGNATURES)
Raised cosine filter.
"""
function rcosfir(β, sps, span = β < 0.68 ? 33-floor(Int, 44β) : 4)
  # default span based on http://www.commsys.isy.liu.se/TSKS04/lectures/3/MichaelZoltowski_SquareRootRaisedCosine.pdf
  # since the span is for rrcosfir, for rcosfir, it is very conservative
  delay = fld(span*sps, 2)
  t = collect(-delay:delay)/sps
  h = Array{Float64}(undef, size(t))
  for i ∈ 1:length(t)
    if abs(t[i]) == 1/(2β)
      h[i] = π/(4sps) * sinc(1/(2β))
    else
      h[i] = sinc(t[i]) * cos(π * β * t[i]) / (1-(2β * t[i])^2) / sps
    end
  end
  h / √sum(h.^2)
end

"""
$(SIGNATURES)
Generates an m-sequence of length `2^m-1` or tap specification `m`.

m-sequences are sequences of `+1/-1` values with near-perfect discrete periodic
auto-correlation properties. All non-zero lag periodic auto-correlations
are -1. The zero-lag autocorrelation is `2^m-1`, where `m` is the shift register
length.

This function currently supports shift register lengths between 2 and 30.

If specification `m` is provided, it should be a list of taps for the shift
register. List of known m-sequence taps can be found in books or
[online](https://users.ece.cmu.edu/~koopman/lfsr/).

# Examples:
```julia-repl
julia> x = mseq(3)                  # generate regular m-sequence
7-element Array{Float64,1}:
  1.0
  1.0
  1.0
 -1.0
  1.0
 -1.0
 -1.0

julia> x = mseq((1,3))              # generate m-sequence with specification (1,3)
7-element Array{Float64,1}:
  1.0
  1.0
  1.0
 -1.0
  1.0
 -1.0
 -1.0
```
"""
function mseq(m, θ=π/2)
  knownspecs = Dict(  # known m-sequences are specified as base 1 taps
       2 => (1,2),          3 => (1,3),          4 => (1,4),          5 => (2,5),
       6 => (1,6),          7 => (1,7),          8 => (1,2,7,8),      9 => (4,9),
      10 => (3,10),        11 => (9,11),        12 => (6,8,11,12),   13 => (9,10,12,13),
      14 => (4,8,13,14),   15 => (14,15),       16 => (4,13,15,16),  17 => (14,17),
      18 => (11,18),       19 => (14,17,18,19), 20 => (17,20),       21 => (19,21),
      22 => (21,22),       23 => (18,23),       24 => (17,22,23,24), 25 => (22,25),
      26 => (20,24,25,26), 27 => (22,25,26,27), 28 => (25,28),       29 => (27,29),
      30 => (7,28,29,30)
  )
  if m ∈ keys(knownspecs)
    spec = collect(knownspecs[m])
  else
    spec = collect(m)
    m = maximum(spec)
  end
  n = 2^m - 1
  reg = ones(UInt8, m)
  x = zeros(Float64, n)
  for j ∈ 1:n
    b = ⊻(reg[spec]...)
    reg = circshift(reg, 1)
    x[j] = 2.0*reg[1] - 1.0
    reg[1] = b
  end
  return x
end

"""
$(SIGNATURES)
Generates an generalized m-sequence of length `2^m-1` or tap specification `m`.

Generalized m-sequences are related to m-sequences but have an additional parameter
`θ`. When `θ = π/2`, generalized m-sequences become normal m-sequences. When
`θ < π/2`, generalized m-sequences contain a DC-component that leads to an exalted
carrier after modulation. When `θ` is `atan(√(2^m-1))`, the m-sequence
is considered to be _period matched_. Period matched m-sequences are complex sequences
with perfect discrete periodic auto-correlation properties, i.e., all non-zero lag
periodic auto-correlations are zero. The zero-lag autocorrelation is `2^m-1`,
where `m` is the shift register length.

This function currently supports shift register lengths between 2 and 30.

# Examples:
```julia-repl
julia> x = gmseq(3)         # generate period matched m-sequence
7-element Array{Complex{Float64},1}:
 0.3535533905932738 + 0.9354143466934853im
 0.3535533905932738 + 0.9354143466934853im
 0.3535533905932738 + 0.9354143466934853im
 0.3535533905932738 - 0.9354143466934853im
 0.3535533905932738 + 0.9354143466934853im
 0.3535533905932738 - 0.9354143466934853im
 0.3535533905932738 - 0.9354143466934853im

julia> x = gmseq(3, π/4)    # generate m-sequence with exalted carrier
7-element Array{Complex{Float64},1}:
 0.7071067811865476 + 0.7071067811865475im
 0.7071067811865476 + 0.7071067811865475im
 0.7071067811865476 + 0.7071067811865475im
 0.7071067811865476 - 0.7071067811865475im
 0.7071067811865476 + 0.7071067811865475im
 0.7071067811865476 - 0.7071067811865475im
 0.7071067811865476 - 0.7071067811865475im
```
"""
function gmseq(m, θ=atan(√(2^maximum(m)-1)))
  x = mseq(m) .+ 0im
  cos(θ) .+ 1im * sin(θ) .* x
end

"""
    hadamard(k)

Generate a Walsh-Hadamard matrix of size `2ᵏ × 2ᵏ`. Each row of the matrix
is orthogonal to all other rows.
"""
function hadamard(k)
  n = 2^k
  [(-1)^count_ones(x&y) for x in 0:n-1, y in 0:n-1]
end

"""
    hadamard(i, k)

Generate a vector with the entries of row `i` of a Walsh-Hadamard matrix of
size `2ᵏ × 2ᵏ`. Rows are numbered from `0` to `2ᵏ-1`, so that `i = 1` is the
first non-trivial (not all ones) Hadamard sequence.
"""
function hadamard(i, k)
  n = 2^k
  0 ≤ i < n || throw(ArgumentError("i must be in the range 0 to 2^k-1"))
  [(-1)^count_ones(x&i) for x in 0:n-1]
end

"""
$(SIGNATURES)
Computes the circular convolution of `x` and `y`. Both vectors must be the same
length.
"""
function circconv(x::AbstractVector, y::AbstractVector)
  length(x) == length(y) || throw(ArgumentError("x and y must be of equal length"))
  ifft(fft(x) .* fft(y))
end

"""
$(SIGNATURES)
Computes the circular correlation of `x` and `y`. Both vectors must be the same
length.
"""
function circcorr(x::AbstractVector, y::AbstractVector=x)
  length(x) == length(y) || throw(ArgumentError("x and y must be of equal length"))
  ifft(conj.(fft(x)) .* fft(y))
end

"""
$(SIGNATURES)
Detects frequency `f` in input signal using the Goertzel algorithm.

The detection metric returned by this function is the complex output
of the Goertzel filter at the end of the input block. Typically, you
would want to compare the magnitude of this output with a threshold to
detect a frequency.

When a block size `n` is specified, the Goertzel algorithm in applied to
blocks of data from the original time series.
"""
function goertzel(x::AbstractVector, f, n; fs=framerate(x))
  signal(map(x1 -> goertzel(x1, f; fs=fs), partition(x, n)), fs / n)
end

function goertzel(x::AbstractVector, f; fs=framerate(x))
  n = length(x)
  m = inHz(f) / (inHz(fs) / n)
  w1 = 0
  w2 = 0
  ω₀ = 2π * m / n
  coeff = 2 * cos(ω₀)
  for j ∈ eachindex(x)
    @inbounds w0 = coeff * w1 - w2 + x[j]
    w2 = w1
    w1 = w0
  end
  w0 = coeff * w1 - w2
  w0 - cis(-ω₀) * w1
end

function goertzel(x::AbstractMatrix, f, n; fs=framerate(x))
  count = cld(size(x, 1), n)
  out = Array{ComplexF64}(undef, (count, nchannels(x)))
  for j ∈ 1:nchannels(x)
    @inbounds out[:, j] = goertzel(x[:, j], f, n; fs=fs)
  end
  signal(out, fs / n)
end

function goertzel(x::AbstractMatrix, f; fs=framerate(x))
  out = Array{ComplexF64}(undef, nchannels(x))
  for j ∈ 1:nchannels(x)
    @inbounds out[j] = goertzel(x[:, j], f; fs=fs)
  end
  out
end

"""
$(SIGNATURES)
Phase-locked loop to track a carrier of frequency approximately `fc` in the input
signal. The loop locks onto carriers within its pull-in range of `fc`, which grows
with `bandwidth`. If `fc` is not specified, the loop starts at 0 Hz.
"""
function pll(x::AbstractVecOrMat, fc=0.0, bandwidth=1e-3; fs=framerate(x))
  fs = inHz(fs)
  fc = inHz(fc)
  β = √bandwidth
  x̄ = samples(x)
  y = Array{ComplexF64}(undef, size(x̄))
  for k ∈ 1:nchannels(x̄)
    xk = @view x̄[:, k]
    rms = √mean(abs2, xk)
    g = rms > 0 ? 1 / rms : 0.0
    ϕ = 0.0
    ω = 0.0
    for j ∈ eachindex(xk)
      y[j, k] = cis(2π * fc * (j - 1) / fs + ϕ)
      Δϕ = angle(xk[j] * conj(y[j, k])) * abs(xk[j]) * g
      ω += bandwidth * Δϕ
      ϕ += β * Δϕ + ω
    end
  end
  signal(y, fs)
end

function sfilt(f::AbstractVector{<:Number}, x::AbstractVector)
  x̄ = samples(x)
  signal(conv(f, x̄)[eachindex(x̄)], framerate(x))
end

sfilt(f, x) = signal(filt(f, samples(x)), framerate(x))
sfilt(b, a, x) = signal(filt(b, a, samples(x)), framerate(x))
sfiltfilt(f, x) = signal(filtfilt(f, samples(x)), framerate(x))
sfiltfilt(b, a, x) = signal(filtfilt(b, a, samples(x)), framerate(x))
sresample(x, rate) = signal(resample(samples(x), rate; dims=1), rate * framerate(x))
sresample(x, rate, coef) = signal(resample(samples(x), rate, coef; dims=1), rate * framerate(x))

"""
    filt(f, x::SampledSignal)
    filt(b, a, x::SampledSignal)

Same as [`filt`](https://docs.juliadsp.org/stable/filters/#DSP.filt),
but retains sampling rate information.
"""
DSP.filt(f::AbstractVector{<:Number}, x::SampledSignal) = sfilt(f, x)
DSP.filt(b::AbstractVector{<:Number}, a::Union{Number,AbstractVector}, x::SampledSignal) = sfilt(b, a, x)

"""
    filtfilt(f, x::SampledSignal)
    filtfilt(b, a, x::SampledSignal)

Same as [`filtfilt`](https://docs.juliadsp.org/stable/filters/#DSP.Filters.filtfilt),
but retains sampling rate information.
"""
DSP.Filters.filtfilt(f::AbstractVector{<:Number}, x::SampledSignal) = sfiltfilt(f, x)
DSP.Filters.filtfilt(b::AbstractVector{<:Number}, a::AbstractVector{<:Number}, x::SampledSignal) = sfiltfilt(b, a, x)

"""
    resample(x::SampledSignal, rate[, coef])

Same as [`resample`](https://docs.juliadsp.org/stable/filters/#DSP.Filters.resample),
but correctly handles sampling rate conversion.
"""
DSP.Filters.resample(x::SampledSignal, rate::Union{Integer,Rational}) = sresample(x, rate)
DSP.Filters.resample(x::SampledSignal, rate::Union{Integer,Rational}, coef::Vector) = sresample(x, rate, coef)

"""
$(SIGNATURES)
Matched filter looking for reference signal `r` in signal `s`.
"""
function mfilter(r, s::AbstractVector)
  issamerate(r, s) || throw(ArgumentError("signals `r` and `s` must have the same sampling rate"))
  if eltype(r) == eltype(s)
    r̄ = samples(r)
    s̄ = samples(s)
  else
    T = promote_type(eltype(r), eltype(s))
    r̄ = convert(AbstractArray{T}, samples(r))
    s̄ = convert(AbstractArray{T}, samples(s))
  end
  f = conj.(reverse(r̄))
  n = length(r) - 1
  sfilt(f, padded(samerateas(s, s̄), (0, n)))[n+1:end]
end

function mfilter(r, s::AbstractMatrix)
  mapreduce(hcat, eachcol(s)) do x
    mfilter(r, x)
  end
end

"""
    istft(Real, X; nfft, noverlap, window)

Compute the inverse short time Fourier transform (ISTFT) of one-sided STFT coefficients `X` which is based
on segments with `nfft` samples with overlap of `noverlap` samples. Refer to `DSP.Periodograms.spectrogram`
for description of the parameters.

For perfect reconstruction, the parameters `nfft`, `noverlap` and `window` in `stft` and
`istft` have to be the same, and the windowing must obey the constraint of "nonzero overlap add" (NOLA).
Implementation based on Zhivomirov 2019 and `istft` in `scipy`.

# Examples:
```julia-repl
julia> x = randn(1024)
1024-element Array{Float64,1}:
 -0.7903319156212055
 -0.564789077302601
  0.8621044972211616
  0.9351928359709288
  ⋮
  2.6158861993992533
  1.2980813993011973
 -0.010592954871694647

julia> X = stft(x, 64, 0)
33×31 Array{Complex{Float64},2}:
  ⋮

julia> x̂ = istft(Real, X; nfft=64, noverlap=0)
1024-element Array{Float64,1}:
 -0.7903319156212054
 -0.5647890773026012
  0.8621044972211612
  0.9351928359709288
  ⋮
  2.6158861993992537
  1.2980813993011973
 -0.010592954871694371
```
"""
function istft(::Type{<:Real}, X::AbstractMatrix{Complex{T}}; nfft::Int, noverlap::Int, window::Union{Function,AbstractVector,Nothing}=nothing) where {T<:AbstractFloat}
  iX = irfft(X, nfft, 1)
  _istft(iX, nfft, noverlap, window)
end

"""
    istft(Complex, X; nfft, noverlap, window)

Compute the inverse short time Fourier transform (ISTFT) of two-sided STFT coefficients `X` which is based
on segments with `nfft` samples with overlap of `noverlap` samples. Refer to `DSP.Periodograms.spectrogram`
for description of the parameters.

For perfect reconstruction, the parameters `nfft`, `noverlap` and `window` in `stft` and
`istft` have to be the same, and the windowing must obey the constraint of "nonzero overlap add" (NOLA).
Implementation based on Zhivomirov 2019 and `istft` in `scipy`.

# Examples:
```julia-repl
julia> x = randn(Complex{Float64}, 1024)
1024-element Array{Complex{Float64},1}:
  -0.5540372432417755 - 0.4286434695080883im
  -0.4759024596520576 - 0.5609424987802376im
                      ⋮
 -0.26493959584225923 - 0.28333817822701457im
  -0.5294529732365809 + 0.7345044619457456im

julia> X = stft(x, 64, 0)
64×16 Array{Complex{Float64},2}:
  ⋮

julia> x̂ = istft(Complex, X; nfft=64, noverlap=0)
1024-element Array{Complex{Float64},1}:
  -0.5540372432417755 - 0.4286434695080884im
 -0.47590245965205774 - 0.5609424987802374im
                      ⋮
  -0.2649395958422591 - 0.28333817822701474im
  -0.5294529732365809 + 0.7345044619457455im
```
"""
function istft(::Type{<:Complex}, X::AbstractMatrix{Complex{T}}; nfft::Int, noverlap::Int, window::Union{Function,AbstractVector,Nothing}=nothing) where {T<:AbstractFloat}
  iX = ifft(X, 1)
  _istft(iX, nfft, noverlap, window)
end

function istft(X::AbstractMatrix{Complex{T}}; nfft::Int, noverlap::Int, window::Union{Function,AbstractVector,Nothing}=nothing) where {T<:AbstractFloat}
  istft(Complex, X; nfft=nfft, noverlap=noverlap, window=window)
end


function _istft(iX::AbstractMatrix{T}, nfft::Int, noverlap::Int, window::Union{Function,AbstractVector,Nothing}=nothing) where {T}
  # H. Zhivomirov, TEM Journal, Vol. 8, No. 1, pp. 56-64, 2019.
  (window === nothing) && (window = rect)
  win, norm2 = Periodograms.compute_window(window, nfft)
  nstep = nfft - noverlap
  nseg = size(iX, 2)
  outputlength = nfft + (nseg-1) * nstep
  iX .*= win
  x = zeros(T, outputlength)
  normw = zeros(eltype(win), outputlength)
  for i = 1:nseg
      @views x[1+(i-1)*nstep:nfft+(i-1)*nstep] .= x[1+(i-1)*nstep:nfft+(i-1)*nstep] .+ iX[:,i]
      @views normw[1+(i-1)*nstep:nfft+(i-1)*nstep] .= normw[1+(i-1)*nstep:nfft+(i-1)*nstep] .+ win .^ 2
  end
  trimlength = nfft % 2 == 0 ? outputlength - nfft : outputlength - nfft + 1
  (sum(@view(normw[1+nfft÷2:end-nfft÷2]) .> 1e-10) != trimlength) && (
    @warn "NOLA condition failed, STFT may not be invertible")
  x .*= nstep/norm2
end

"""
$(SIGNATURES)
Spectral whitening of input signal `x` in the frequency domain. The parameters `nfft`,
`noverlap` and `window` are required for the computation of STFT coefficients of `x`.
Refer to `DSP.Periodograms.spectrogram` for description of the parameters. `γ` is a
scaling or degree-of-flattening factor. The algorithm is based on Lee 1986.
"""
function whiten(x::AbstractVector; nfft::Int, noverlap::Int, window::Union{Function,AbstractVector,Nothing}=nothing, γ=1)
  # M. W. Lee, Open-File Report 86-108, 1986.
  xstft = stft(x, nfft, noverlap; window=window)
  mag = abs.(xstft)
  logmag = log.(mag .+ eps(eltype(mag)))
  logmag .-= γ * mean(logmag; dims=2)
  outputtype = isreal(x) ? Real : Complex
  istft(outputtype, exp.(logmag) .* exp.(im .* angle.(xstft)); nfft=nfft, noverlap=noverlap, window=window)
end

"""
    delay!(x, v)

Delay signal `x` by `v` units. The delay may be specified in samples or in time
units (e.g. 2.3𝓈). Delays may be negative, and may be a fractional number of
samples. For a multichannel signal, every channel is delayed. The length of the
returned signal is the same as the original, with any part of the signal that is
shifted beyond the original length discarded, and zeros shifted in.

Fractional delays use the same band-limited interpolation as evaluating a signal
at a time (see [`SampledSignal`](@ref)), so a signal `s` delayed by `τ` seconds
has the value `s(t - τ)` at every sample time `t`. Frequencies above about 90% of
the Nyquist frequency are attenuated.

# Examples:
```julia-repl
julia> x = signal([1.0, 2.0, 3.0], 10.0)
julia> delay!(x, 1)
SampledSignal @ 10.0 Hz, 3-element Vector{Float64}:
 0.0
 1.0
 2.0

julia> x = signal([1.0, 2.0, 3.0], 10.0)
julia> delay!(x, -1)
SampledSignal @ 10.0 Hz, 3-element Vector{Float64}:
 1.0
 2.0
 0.0

julia> x = signal([1.0, 2.0, 3.0, 2.0, 1.0], 10.0)
julia> delay!(x, 1.2)
SampledSignal @ 10.0 Hz, 5-element Vector{Float64}:
 -0.10392666815161294
  0.8034767135536652
  1.7503401027230303
  2.9432140000798874
  2.274536909816327

julia> x = signal([1.0, 2.0, 3.0, 2.0, 1.0], 10.0)
julia> delay!(x, -0.11𝓈)
SampledSignal @ 10.0 Hz, 5-element Vector{Float64}:
  2.1354604073929195
  2.985636168304411
  1.8710195275583208
  0.9043683839765194
 -0.06009907434381305
```
"""
function delay!(x, n::Integer)
  for j ∈ axes(x, 2)
    shift!(view(x, :, j), n)
  end
  x
end

function delay!(x, v::Real)
  vi = floor(Int, v)
  v == vi && return delay!(x, vi)
  if abs(v) > size(x, 1) + SINC_HALFWIDTH
    x .= 0
    return x
  end
  w = sincweights(real(float(eltype(x))), 1 - (v - vi))
  for j ∈ axes(x, 2)
    c = view(x, :, j)
    p = padded(copy(samples(c)))
    for i ∈ eachindex(c)
      c[i] = sincsum(p, i - vi - 1, w, 1)
    end
  end
  x
end

delay!(x, t::Units.Unitful.Time) = delay!(x, inseconds(t) * framerate(x))

"""
    delay(x, v)

Create a delayed version of signal `x` with `v` units of delay.
See [`delay!`](@ref) for details.
"""
delay(x, v) = delay!(copy(x), v)

function shift!(x, n)
  if n > length(x) || n < -length(x)
    x .= 0
  elseif n > 0
    @inbounds for i ∈ lastindex(x):-1:firstindex(x)+n
      x[i] = x[i-n]
    end
    x[begin:begin+n-1] .= 0
  elseif n < 0
    @inbounds for i ∈ firstindex(x):lastindex(x)+n
      x[i] = x[i-n]
    end
    x[end+n+1:end] .= 0
  end
  x
end

"""
$(SIGNATURES)
Finds up to `n` strongest copies of reference signal `r` in signal `s`. The reference
signal `r` should have a delta-like autocorrelation for this function to work
well.

If the keyword parameter `finetune` is set to 0, approximate arrival times are
computed based on a matched filter. If it is set to a positive value, an
iterative optimization is performed to find more accruate arrival times within
roughly `finetune` samples of the coarse arrival times.

The keyword parameter `mingap` controls the minimum number of samples between
two arrivals.

Returns named tuple `(time=t, amplitude=a, mfo=m)` where `t` is a vector of arrival
times, `a` is a vector of complex amplitudes of the arrivals, and `m` is the complex
matched filter output (if keyword parameter `mfo` is set to `true`). The arrivals
are sorted in ascending order of arrival times.

# Examples:
```julia-repl
julia> x = chirp(1000, 5000, 0.1, 40960; window=(tukey, 0.05))
julia> x4 = resample(x, 4)
julia> y4 = samerateas(x4, zeros(32768))
julia> y4[128:127+length(x4)] = real(x4)          # time 0.000775𝓈, index 32.75
julia> y4[254:253+length(x4)] += -0.8 * real(x4)  # time 0.001544𝓈, index 64.25
julia> y4[513:512+length(x4)] += 0.6 * real(x4)   # time 0.003125𝓈, index 129.0
julia> y = resample(y4, 1//4)
julia> y .+= 0.1 * randn(length(y))
julia> findsignal(x, y, 3)
(time = [0.000781, 0.001538, 0.003125], amplitude = [...], mfo=[])
julia> findsignal(x, y, 3; mfo=true)
(time = [0.000775, 0.001545, 0.003124], amplitude = [...], mfo=[...])
```
"""
# golden-section minimization of a unimodal function over [a, b]
function _minimize(f, a, b; tol=1e-4)
  φ = (√5 - 1) / 2
  c = b - φ * (b - a)
  d = a + φ * (b - a)
  fc = f(c)
  fd = f(d)
  while b - a > tol
    if fc < fd
      b, d, fd = d, c, fc
      c = b - φ * (b - a)
      fc = f(c)
    else
      a, c, fc = c, d, fd
      d = a + φ * (b - a)
      fd = f(d)
    end
  end
  (a + b) / 2
end

function findsignal(r, s, n=1; prominence=0.0, finetune=2, mingap=1, mfo=false)
  # coarse arrival time estimation
  r = analytic(r)
  s = analytic(s)
  m = mfilter(r, s)
  m ./= (std(r) * length(r))
  T = eltype(m)
  m̄ = abs.(samples(m))
  p = argmaxima(m̄, mingap)
  prominence > 0 && peakproms!(p, m̄; min=prominence*maximum(m̄))
  if length(p) > length(s)/10
    return (time=Float64[], amplitude=T[], mfo=mfo ? m : empty(m))
  end
  h = m̄[p]
  ndx = partialsortperm(h, 1:min(n,length(h)); rev=true)
  p = p[ndx]
  if finetune == 0
    t = time(Float64.(p), s)
    ndx = sortperm(t)
    return (time=t[ndx], amplitude=samples(m[p[ndx]]), mfo=mfo ? m : empty(m))
  end
  # iterative fine arrival time estimation
  i::Int = minimum(p)
  N::Int = maximum(p) - i + length(r) + 2 * finetune
  N = nextfastfft(N)
  i = max(1, i - finetune)
  X = fft(vcat(samples(r), zeros(N-length(r))))
  S = fft(vcat(samples(s[i:min(i+N-1,end)]), zeros(max(i+N-1-length(s),0))))
  ii, a = let P=length(p), f=fftfreq(N)
    # variable projection: for fixed fractional delays, the optimal amplitudes
    # are given by a linear least-squares solve, so we only need to search over
    # the delays (each bracketed within ±finetune samples of its coarse peak)
    ii = Float64.(p .- i)
    B = Matrix{ComplexF64}(undef, N, P)
    for k ∈ 1:P
      @views B[:,k] .= X .* cis.(-2π .* ii[k] .* f)
    end
    cost = () -> (a = B \ S; sum(abs2, B * a .- S))
    for _ ∈ 1:3
      for k ∈ 1:P
        ii[k] = _minimize(ii[k] - finetune, ii[k] + finetune; tol=1e-3) do v
          @views B[:,k] .= X .* cis.(-2π .* v .* f)
          cost()
        end
        @views B[:,k] .= X .* cis.(-2π .* ii[k] .* f)
      end
    end
    (ii, B \ S)
  end
  pp = ii .+ i
  t = time(pp, s)
  ndx = sortperm(t)
  (time=t[ndx], amplitude=a[ndx], mfo=mfo ? m : empty(m))
end

"""
$(SIGNATURES)
Compose a signal from a reference signal and a list of arrival times and amplitudes.

# Examples:
```julia-repl
julia> x = cw(10kHz, 0.01, 44.1kHz)
julia> y1 = compose(x, [0.01, 0.03, 0.04], [1.0, 0.8, 0.6]; duration=0.05)
julia> y2 = compose(real(x), [10ms, 30ms, 40ms], [1.0, 0.8, 0.6]; duration=50ms)
```
"""
function compose(r, t, a; duration=duration(r)+maximum(t), fs=framerate(r))
  r isa SampledSignal && framerate(r) != fs && throw(ArgumentError("Reference signal must be sampled at fs"))
  length(t) == length(a) || @warn "Mismatch in number of time and amplitude entries, using minimum of the two..."
  r̃ = analytic(r)
  x = zeros(eltype(r̃), ceil(Int, inseconds(duration) * inHz(fs)))
  x1 = copy(x)
  for (t1, a1) ∈ zip(t, a)
    x1 .= 0
    copyto!(x1, 1, samples(r̃), 1, min(length(x1), nframes(r̃)))
    delay!(x1, t1 * inHz(fs))
    x1 .*= a1
    x .+= x1
  end
  signal(isanalytic(r) ? x : √2 * real(x), fs)
end

"""
$(SIGNATURES)
Decompose a signal as a sum of reference signal with some arrival times and amplitudes.
If the number of signals `n` is non-zero, the function will limit the decomposition
to the strongest `n` signals.

The resolution for arrival time is one sample. The `threshold` parameter controls
the stopping criterion for the decomposition. If the relative change in the residual
energy is less than `threshold`, the decomposition stops. The `refine` parameter
controls whether the amplitudes are refined using an optimization algorithm. Without
refinement, the amplitudes are simply the matched filter outputs at the arrival times.

# Examples:
```julia-repl
julia> x = compose(mseq(12), [0.1, 0.2], [1.0, 0.7]; duration=1.0, fs=8000)
julia> x += 0.1 * randn(size(x))
julia> decompose(mseq(12), x)
(time=[0.1, 0.2], amplitude=[1.00201, 0.70061], index=[801, 1601])
```
"""
function decompose(r::AbstractVector, x::AbstractVector, n=0; threshold=0.01, refine=true)
  # use OMP algorithm for sparse signal decomposition
  fs = framerate(x)
  Λ = Int[]
  a = Array{eltype(x)}(undef, 0)
  xᵣ = x
  while n == 0 || length(Λ) < n
    mfo = mfilter(r, xᵣ) / energy(r)
    absmfo = abs.(mfo)
    absmfo[Λ] .= 0
    λ = argmax(absmfo)
    absmfo[λ] == 0 && break
    push!(Λ, λ)
    push!(a, mfo[λ])
    if refine
      # amplitudes with fixed arrival times form a linear least-squares problem
      A = stack(λ -> samples(compose(r, [(λ - 1) / fs], [one(eltype(x))]; duration=duration(x), fs)), Λ)
      if eltype(A) <: Real && eltype(x) <: Complex
        # real reference with complex amplitudes: compose output √2⋅Re(Σ aₖr̃ₖ) is
        # only R-linear in a, so solve for real & imaginary parts separately
        B = stack(λ -> samples(compose(r, [(λ - 1) / fs], [im * one(eltype(x))]; duration=duration(x), fs)), Λ)
        v = [A B] \ real(samples(x))
        next_a = complex.(v[1:length(Λ)], v[length(Λ)+1:end])
      else
        next_a = A \ samples(x)
      end
    else
      next_a = a
    end
    x̂ = compose(r, (Λ .- 1) ./ fs, next_a; duration=duration(x), fs)
    prev_rms = rms(xᵣ[λ:min(λ+length(r)-1, end)])
    xᵣ = x - x̂
    if (prev_rms - rms(xᵣ[λ:min(λ+length(r)-1, end)])) / prev_rms < threshold
      # remove last entry, since it does not meet threshold
      resize!(Λ, length(Λ)-1)
      resize!(a, length(a)-1)
      break
    end
    a = next_a
  end
  ndx = sortperm(Λ)
  (time=(Λ[ndx] .- 1) ./ fs, amplitude=a[ndx], index=Λ[ndx])
end

"""
    dzt(x, L, K)
    dzt(x, L)

Compute the discrete Zak transform (DZT) of signal `x` with `L` delay bins and `K`
Doppler bins. The length of signal `x` must be equal to `LK`. If `K` is not specified,
it is assumed to be the length of `x` divided by `L`. Returns a `K × L` complex matrix.

If the frame rate of `x` if `fs` Sa/s, the delay bins are spaced at `1/fs` seconds
and the Doppler bins are spaced at `fs/LK` Hz. The DZT is scaled such
that the energy of the signal is preserved, i.e., `sum(abs2, x) ≈ sum(abs2, X)`.

For efficient computation of the DZT, `K` should product of small primes.

# Examples:
```julia-repl
julia> x = randn(ComplexF64, 4096)
4096-element Vector{ComplexF64}:
  :
julia> X = dzt(x, 64)
64×64 Matrix{ComplexF64}:
  :
```
"""
function dzt(x::AbstractVector, L::Int, K::Int)
  length(x) == L * K || throw(ArgumentError("Length of x must be L * K"))
  X = complex.(collect(transpose(reshape(x, L, K)))) ./ sqrt(K)
  fft!(X, 1)
end

function dzt(x::AbstractVector, L::Int)
  length(x) % L == 0 || throw(ArgumentError("Length of x must be a multiple of L"))
  dzt(x, L, length(x) ÷ L)
end

"""
$(SIGNATURES)
Compute the inverse discrete Zak transform (DZT) of 2D `K × L` complex signal `X`
with `L` delay bins and `K` Doppler bins. The length of the returned signal is `LK`.

See [`dzt`](@ref) for more details.

# Examples:
```julia-repl
julia> x = randn(ComplexF64, 4096)
4096-element Vector{ComplexF64}:
  :
julia> X = dzt(x, 64)
64×64 Matrix{ComplexF64}:
  :
julia> idzt(X) ≈ x
true
```
"""
function idzt(X::AbstractMatrix)
  X = ifft(complex.(X), 1)
  X .*= sqrt(size(X, 1))
  collect(vec(transpose(X)))
end
