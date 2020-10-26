import DSP: filt, filtfilt, resample

export fir, removedc, removedc!, demon
export upconvert, downconvert, rrcosfir, rcosfir
export mseq, gmseq, circconv, goertzel, pll
export sfilt, sfiltfilt, sresample, mfilter
export istft, spectralwhitening

"""
$(SIGNATURES)
Designs a `n`-tap FIR filter with a passband from `f1` to `f2` using the
specified `method`. If frame rate `fs` is not specified, `f1` and `f2` are given
in normalized units (1.0 being Nyquist). If `f1` is 0, the designed filter is
a lowpass filter, and if `f2` is `nothing` then it is a highpass filter.

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
```
"""
function fir(n, f1, f2=nothing; fs=2.0, method=FIRWindow(hanning(n)))
  fs = inHz(fs)
  if f1 == 0
    f = Lowpass(inHz(f2); fs=fs)
  elseif f2 == nothing || inHz(f2) == fs/2
    f = Highpass(inHz(f1); fs=fs)
  else
    f = Bandpass(inHz(f1), inHz(f2); fs=fs)
  end
  return digitalfilter(f, method)
end

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
  pad = cld(length(pulseshape), 2*sps) - 1
  s = vcat(zeros(pad), analytic(s), zeros(pad))
  s = signal(resample(s, sps, pulseshape), sps*fs)
  √2 * real.(s .* cis.(2π * inHz(fc) * domain(s)))
end

function upconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  out = Any[]
  for j in 1:nchannels(s)
    push!(out, upconvert(s[:,j], sps, fc, pulseshape; fs=fs))
  end
  hcat(out...)
end

"""
$(SIGNATURES)
Converts passband signal centered around carrier frequency `fc` to baseband,
and downsamples it by a factor of `sps`. If the `pulseshape` is specified to
be `nothing`, downsampling is performed without filtering.
"""
function downconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  s = signal(analytic(s), fs)
  s .*= cis.(-2π * inHz(fc) * domain(s))
  sps == 1 && return signal(s, fs)
  pulseshape == nothing && return signal(s[1:sps:end,:], fs/sps)
  signal(resample(s, 1//sps, pulseshape), fs/sps)
end

function downconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  out = Any[]
  for j in 1:nchannels(s)
    push!(out, downconvert(s[:,j], sps, fc, pulseshape; fs=fs))
  end
  hcat(out...)
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
$(SIGNATURES)
Computes the circular convolution of `x` and `y`. Both vectors must be the same
length.
"""
function circconv(x::AbstractVector, y::AbstractVector=x)
  if length(x) != length(y)
    throw(ArgumentError("x and y must be of equal length"))
  end
  n = length(x)
  z = similar(x)
  for j ∈ 1:n
    z[j] = circshift(x, j-1)' * y
  end
  return z
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
  signal(map(x1 -> goertzel(x1, f; fs=fs), partition(x, n)), fs/n)
end

function goertzel(x::AbstractVector, f; fs=framerate(x))
  n = length(x)
  m = inHz(f)/(inHz(fs)/n)
  w1 = 0
  w2 = 0
  for j ∈ 1:n
    w0 = 2 * cos(2π * m/n) * w1 - w2 + x[j]
    w2 = w1
    w1 = w0
  end
  w0 = 2 * cos(2π * m/n) * w1 - w2
  w0 - cis(-2π * m/n) * w1
end

function goertzel(x::AbstractMatrix, f, n; fs=framerate(x))
  count = cld(size(x,1), n)
  out = Array{ComplexF64}(undef, (count, nchannels(x)))
  for j ∈ 1:nchannels(x)
    out[:,j] = goertzel(x[:,j], f, n; fs=fs)
  end
  signal(out, fs/n)
end

function goertzel(x::AbstractMatrix, f; fs=framerate(x))
  out = Array{ComplexF64}(undef, nchannels(x))
  for j ∈ 1:nchannels(x)
    out[j] = goertzel(x[:,j], f; fs=fs)
  end
  out
end

"""
$(SIGNATURES)
Phased-lock loop to track dominant carrier frequency in input signal.
"""
function pll(x::AbstractVecOrMat, bandwidth=1e-3; fs=framerate(x))
  β = √bandwidth
  n = nchannels(x)
  ϕ = zeros(1,n)
  ω = zeros(1,n)
  y = similar(x, ComplexF64)
  for j ∈ 1:nframes(x)
    y[j,:] = cis.(ϕ)
    Δϕ = angle.(x[j,:] .* conj.(y[j,:])) .* abs.(x[j,:])
    ω .+= bandwidth * Δϕ
    ϕ .+= β*Δϕ .+ ω
  end
  signal(y, fs)
end

"""
$(SIGNATURES)
Same as [`filt`](https://docs.juliadsp.org/stable/filters/#DSP.filt),
but retains sampling rate information.
"""
sfilt(f, x, args...) = signal(filt(f, samples(x), args...), framerate(x))
sfilt(b, a, x, args...) = signal(filt(b, a, samples(x), args...), framerate(x))

"""
$(SIGNATURES)
Same as [`filtfilt`](https://docs.juliadsp.org/stable/filters/#DSP.Filters.filtfilt),
but retains sampling rate information.
"""
sfiltfilt(coef, x) = signal(filtfilt(coef, samples(x)), framerate(x))

"""
$(SIGNATURES)
Same as [`resample`](https://docs.juliadsp.org/stable/filters/#DSP.Filters.resample),
but correctly handles sampling rate conversion.
"""
sresample(x, rate, args...) = signal(resample(samples(x), rate, args...), rate * framerate(x))

"""
$(SIGNATURES)
Matched filter looking for reference signal `r` in signal `s`.
"""
function mfilter(r, s)
  f = conj.(reverse(samples(r)))
  n = length(r) - 1
  sfilt(f, padded(s, (0, n)))[n+1:end]
end

"""
$(SIGNATURES)
Compute the inverse short time Fourier transform (ISTFT) of STFT coefficients `X` which is based
on segments with `n` samples with overlap of `noverlap` samples. Refer to `DSP.Periodograms.spectrogram` 
for description of optional keyword arguments.

For perfect reconstruction, the parameters `n`, `noverlap` and `window` in `stft` and 
`istft` have to be the same, and the windowing must obey the constraint of "nonzero overlap add" (NOLA).  
Implementation based on Zhivomirov 2019 and `istft` in `scipy`. 

H. Zhivomirov. On the Development of STFT-analysis and ISTFT-synthesis Routines 
and their Practical Implementation. TEM Journal, ISSN: 2217-8309, DOI: 10.18421/TEM81-07, 
Vol. 8, No. 1, pp. 56-64, Feb. 2019. 
(http://www.temjournal.com/content/81/TEMJournalFebruary2019_56_64.pdf)

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
julia> x̂ = istft(X, 64, 0)
1024-element Array{Float64,1}:
 -0.7903319156212054
 -0.5647890773026012
  0.8621044972211612
  0.9351928359709288
  ⋮
  2.6158861993992537
  1.2980813993011973
 -0.010592954871694371
"""
function istft(X::AbstractMatrix{Complex{T}}, 
               n::Int, 
               noverlap::Int; 
               onesided::Bool=true, 
               window::Union{Function,AbstractVector,Nothing}=nothing) where {T<:AbstractFloat}
  (window === nothing) && (window = rect)
  win, norm2 = Periodograms.compute_window(window, n)

  nstep = n - noverlap
  nseg = size(X, 2)
  outputlength = n + (nseg-1) * nstep

  iX = onesided ? irfft(X, n, 1) : ifft(X, 1)
  iX .*= win
  x = zeros(eltype(iX), outputlength)
  normw = zeros(T, outputlength)
  for i = 1:nseg
      x[1+(i-1)*nstep:n+(i-1)*nstep] += iX[:,i]
      normw[1+(i-1)*nstep:n+(i-1)*nstep] += win .^ 2
  end

  (sum(normw[n÷2:end-n÷2] .> 1e-10) != length(normw[n÷2:end-n÷2])) && (
      @warn "NOLA condition failed, STFT may not be invertible")
  x .*= nstep/norm2
end

"""
$(SIGNATURES)
Spectral whitening of input signal `x`. The parameters `n`, `noverlap` and 
`window` are required for the computation of STFT coefficients of `x`. Refer to 
`DSP.Periodograms.spectrogram` for description of the parameters. `γ` is a scaling 
or degree-of-flattening factor. The algorithm is based on Lee 1986.

M. W. Lee, Spectral whitening in the frequency domain, 1986.
"""
function whiten(x::AbstractVector, 
                n::Int, 
                noverlap::Int; 
                window::Union{Function,AbstractVector,Nothing}=nothing,
                γ=1)
  xstft = stft(x, n, noverlap; window=window)
  mag = abs.(xstft)
  logmag = log.(mag .+ eps(eltype(mag)))
  logmag .-= γ * mean(logmag; dims=2)
  istft(exp.(logmag) .* exp.(im .* angle.(xstft)), n, noverlap; window=window) 
end