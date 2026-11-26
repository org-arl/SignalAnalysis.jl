using FFTW: AbstractFFTs

export spectrum, tobin

struct SpectrumInfo
  fs::Float32
  n::Int
  onesided::Bool
  k0::Int
end

Base.keys(x::SpectrumInfo) = fieldnames(typeof(x))

const SampledSpectrum = MetaArray{<:Any,SpectrumInfo,T} where T

struct PeriodicView{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  data::A
  info::SpectrumInfo
end

const PaddedSpectrum = MetaArray{<:PeriodicView,SpectrumInfo,T} where T
const NumericSpectrum = Union{MetaArray{<:Any,SpectrumInfo,<:Number,1},MetaArray{<:Any,SpectrumInfo,<:Number,2}}

withinfo(X, m::SpectrumInfo) = MetaArray(m, X)

spectruminfo(X::SampledSpectrum) = getfield(X, :meta)

firstbin(X) = spectruminfo(X).k0 + firstindex(X, 1) - 1
lastbin(X) = spectruminfo(X).k0 + lastindex(X, 1) - 1

"""
    spectrum(H, fs; onesided=false)
    spectrum(H, fs; onesided=true, n)

Creates a spectrum from FFT coefficients `H` of a signal sampled at rate `fs`.
Each column of `H` is a channel. A two-sided spectrum holds all `n = size(H, 1)`
bins of an `n`-point FFT. A one-sided spectrum holds the `n÷2+1` non-negative
frequency bins; `n` is only needed if it is odd.

# Examples:
```julia-repl
julia> h = ifft(spectrum(H, 8kHz))   # impulse response from frequency response H
```
"""
function spectrum(H::AbstractArray, fs; onesided=false, n=nothing)
  m = size(H, 1)
  if onesided
    n === nothing && (n = max(2 * (m - 1), 1))
    n ÷ 2 + 1 == m || throw(ArgumentError("a one-sided spectrum of an $n-point FFT has $(n ÷ 2 + 1) bins, not $m"))
  else
    n === nothing || throw(ArgumentError("n is only needed for one-sided spectra"))
    n = m
  end
  withinfo(H, SpectrumInfo(inHz(fs), n, onesided, 0))
end

fftspectrum(X, x::SampledSignal; onesided=false) =
  withinfo(X, SpectrumInfo(framerate(x), nframes(x), onesided, 0))

SignalBase.framerate(X::SampledSpectrum) = X.fs

samples(X::SampledSpectrum) = samples(getfield(X, :data))

function Base.show(io::IO, mime::MIME"text/plain", X::SampledSpectrum)
  m = spectruminfo(X)
  print(io, "SampledSpectrum of ", m.n, "-point signal @ ", m.fs, " Hz, bins ", firstbin(X), ":", lastbin(X))
  m.onesided && print(io, ", one-sided")
  X isa PaddedSpectrum && print(io, ", padded")
  print(io, ", ")
  show(io, mime, samples(X))
end

### periodic view of padded spectra

Base.size(p::PeriodicView) = size(p.data)
Base.axes(p::PeriodicView) = axes(p.data)

Base.checkbounds(::Type{Bool}, p::PeriodicView, i::PaddedIndex, j::PaddedIndex...) =
  Base.checkbounds_indices(Bool, Base.tail(axes(p.data)), j)

function binvalue(p::PeriodicView{T}, k, j...) where T
  m = p.info
  b1 = m.k0 + firstindex(p.data, 1) - 1
  b1 ≤ k ≤ b1 + size(p.data, 1) - 1 ? p.data[firstindex(p.data, 1) + k - b1, j...] : zero(T)
end

function periodicvalue(p::PeriodicView, i, j...)
  m = p.info
  k = m.k0 + i - 1
  if m.onesided
    kk = mod(k, m.n)
    kk ≤ m.n ÷ 2 ? binvalue(p, kk, j...) : conj(binvalue(p, m.n - kk, j...))
  else
    b1 = m.k0 + firstindex(p.data, 1) - 1
    binvalue(p, b1 + mod(k - b1, m.n), j...)
  end
end

Base.getindex(p::PeriodicView{T,1}, i::Int) where T = periodicvalue(p, i)
Base.getindex(p::PeriodicView{T,2}, i::Int, j::Int) where T = periodicvalue(p, i, j)

samples(p::PeriodicView) = samples(p.data)

### padding

"""
    padded(X::SampledSpectrum)
    padded(X::SampledSpectrum, padding)

Pads spectrum `X` with zeros. Without `padding`, the spectrum is padded
indefinitely: the result is a view with the same indices as `X`, but any index is
valid, giving zeros within one period of the FFT grid, and repeating periodically
beyond it. One-sided spectra are extended by conjugate symmetry, as for a real
signal.

With `padding` (a number, or a tuple of left and right amounts), the result is a
new spectrum indexed from 1, extended by that many zero bins on each side, within
one period of its FFT grid. Its frequencies stay correct, since the spectrum keeps
track of its bins.

A padded spectrum can be inverted with `ifft` or `irfft`, which select one period
of it. This allows a band of a spectrum to be inverted, e.g. `ifft(padded(X[a:b]))`.
"""
nokeywords(kwargs) = isempty(kwargs) || throw(ArgumentError("spectra are padded with zeros, without keyword options"))

padded(X::PaddedSpectrum; kwargs...) = (nokeywords(kwargs); X)

function padded(X::SampledSpectrum; kwargs...)
  nokeywords(kwargs)
  withinfo(PeriodicView(getfield(X, :data), spectruminfo(X)), spectruminfo(X))
end

function padded(X::SampledSpectrum, padding; kwargs...)
  nokeywords(kwargs)
  X isa PaddedSpectrum && throw(ArgumentError("spectrum is already padded indefinitely"))
  left, right = length(padding) == 1 ? (padding, padding) : (padding[1], padding[2])
  m = spectruminfo(X)
  if m.onesided
    firstbin(X) - left ≥ 0 && lastbin(X) + right ≤ m.n ÷ 2 ||
      throw(ArgumentError("one-sided spectrum can only be padded within bins 0:$(m.n ÷ 2)"))
  else
    size(X, 1) + left + right ≤ m.n ||
      throw(ArgumentError("spectrum can only be padded within one period of $(m.n) bins"))
  end
  withinfo(samples(padded(samples(X), (left, right))), SpectrumInfo(m.fs, m.n, m.onesided, firstbin(X) - left))
end

### indexing

toindex(X, i) = i

function toindex(X, t::Tuple{Units.Unitful.Frequency,Units.Unitful.Frequency})
  r = tobin(t, X)
  X isa PaddedSpectrum || Base.checkindex(Bool, axes(X, 1), r) ||
    throw(ArgumentError("frequency band $(t[1]) to $(t[2]) is outside the spectrum; use fftshift() or padded() for negative frequencies"))
  r
end

sliceinfo(X, i) = ndims(X) == 1 ? sliceinfo(spectruminfo(X), i) : nothing
sliceinfo(X, i, j...) = sliceinfo(spectruminfo(X), i)
sliceinfo(m::SpectrumInfo, i) = nothing
sliceinfo(m::SpectrumInfo, ::Colon) = m
sliceinfo(m::SpectrumInfo, i::AbstractUnitRange{<:Integer}) = SpectrumInfo(m.fs, m.n, m.onesided, m.k0 + first(i) - 1)

function sliceinfo(m::SpectrumInfo, i::AbstractRange{<:Integer})
  s = step(i)
  s == 1 && return SpectrumInfo(m.fs, m.n, m.onesided, m.k0 + first(i) - 1)
  k = m.k0 + first(i) - 1
  s > 0 && mod(m.n, s) == 0 && mod(k, s) == 0 ||
    throw(ArgumentError("stride $s does not align with the FFT grid; use samples() for arbitrary bin selections"))
  SpectrumInfo(m.fs, m.n ÷ s, m.onesided, k ÷ s)
end

function spectrumindex(f, X, I...)
  J = map(i -> toindex(X, i), I)
  data = getfield(X, :data)
  m = sliceinfo(X, J...)
  m === nothing ? f(data, J...) : withinfo(f(data, J...), m)
end

Base.getindex(X::SampledSpectrum, i...) = spectrumindex(getindex, X, i...)
Base.getindex(X::SampledSpectrum, i::Int...) = getindex(getfield(X, :data), i...)
Base.view(X::SampledSpectrum, i...) = spectrumindex(view, X, i...)
Base.view(X::SampledSpectrum, i::Union{Real,AbstractArray}...) = spectrumindex(view, X, i...)
Base.view(X::SampledSpectrum, i::Union{Real,Base.AbstractCartesianIndex,AbstractArray}...) = spectrumindex(view, X, i...)

Base.setindex!(X::SampledSpectrum, v, t::Tuple{Units.Unitful.Frequency,Units.Unitful.Frequency}, j...) =
  setindex!(X, v, tobin(t, X), j...)

### concatenation

function samegrid(X, Y)
  a, b = spectruminfo(X), spectruminfo(Y)
  a.fs == b.fs && a.n == b.n && a.onesided == b.onesided
end

function Base.vcat(X::NumericSpectrum, Xs::AbstractVecOrMat{<:Number}...)
  all(x -> x isa SampledSpectrum, Xs) || return vcat(samples(X), map(samples, Xs)...)
  prev = X
  for x ∈ Xs
    samegrid(X, x) && firstbin(x) == lastbin(prev) + 1 ||
      throw(ArgumentError("spectra must be contiguous segments of the same FFT grid to be concatenated"))
    prev = x
  end
  m = spectruminfo(X)
  m.onesided || lastbin(prev) - firstbin(X) < m.n ||
    throw(ArgumentError("concatenated spectrum exceeds one period of $(m.n) bins"))
  withinfo(vcat(samples(X), map(samples, Xs)...), SpectrumInfo(m.fs, m.n, m.onesided, firstbin(X)))
end

function Base.hcat(X::NumericSpectrum, Xs::AbstractVecOrMat{<:Number}...)
  all(x -> x isa SampledSpectrum, Xs) || return hcat(samples(X), map(samples, Xs)...)
  all(x -> samegrid(X, x) && firstbin(x) == firstbin(X) && lastbin(x) == lastbin(X), Xs) ||
    throw(ArgumentError("spectra must have identical FFT grids to be concatenated as channels"))
  m = spectruminfo(X)
  withinfo(hcat(samples(X), map(samples, Xs)...), SpectrumInfo(m.fs, m.n, m.onesided, firstbin(X)))
end

### frequencies

"""
    domain(X::SampledSpectrum)

Returns the frequencies (in Hz) of the bins of spectrum `X`, as a range. The
frequencies are not wrapped, so an unshifted two-sided spectrum spans `0` to
`fs`, and a shifted one spans `-fs/2` to `fs/2`. Use [`fftfreq`](@ref) for
frequencies wrapped to the interval `[-fs/2, fs/2)`.
"""
domain(X::SampledSpectrum) = (firstbin(X):lastbin(X)) .* (framerate(X) / spectruminfo(X).n)

"""
    fftfreq(X::SampledSpectrum)

Returns the frequencies (in Hz) of the bins of spectrum `X`, wrapped to the
interval `[-fs/2, fs/2)` for two-sided spectra.
"""
function AbstractFFTs.fftfreq(X::SampledSpectrum)
  f = collect(domain(X))
  spectruminfo(X).onesided && return f
  fs = framerate(X)
  mod.(f .+ fs / 2, fs) .- fs / 2
end

"""
    fftfreq(x::SampledSignal)
    fftfreq(n, x::SampledSignal)

Frequencies of the bins of an `n`-point FFT of signal `x`. If `n` is not
specified, the number of frames in `x` is used.
"""
AbstractFFTs.fftfreq(x::SampledSignal) = fftfreq(nframes(x), framerate(x))
AbstractFFTs.fftfreq(n::Integer, x::SampledSignal) = fftfreq(n, framerate(x))

"""
    tobin(f, X::SampledSpectrum)
    tobin(f, fs, n)

Converts frequency `f` to the index of the nearest bin in spectrum `X`, or in an
unshifted `n`-point FFT at sampling rate `fs`. Frequencies without units are in Hz.
Like [`toframe`](@ref), `f` may be a vector, a range or a tuple giving a frequency
band. Indices are not wrapped or bounds checked, so negative frequencies in an
unshifted spectrum map to indices that are only valid in a padded spectrum.

# Examples:
```julia-repl
julia> X = fft(signal(randn(1000), 8kHz));
julia> tobin(1kHz, X)
126

julia> tobin((1kHz, 2kHz), X)
126:251
```
"""
tobin(f, X::SampledSpectrum) = round.(Int, inHz.(f) .* (spectruminfo(X).n / Float64(framerate(X)))) .- spectruminfo(X).k0 .+ 1
tobin(f::Tuple{Any,Any}, X::SampledSpectrum) = UnitRange(tobin(f[1], X), tobin(f[2], X))
tobin(f, fs, n) = round.(Int, inHz.(f) .* (n / inHz(fs))) .+ 1
tobin(f::Tuple{Any,Any}, fs, n) = UnitRange(tobin(f[1], fs, n), tobin(f[2], fs, n))

"""
    freq(i, X::SampledSpectrum)

Gets the frequency (in Hz) of the `i`th bin of spectrum `X`. The index `i` can be
a range or an array of indices, and may lie outside `X` for padded spectra.
"""
DSP.freq(i::Real, X::SampledSpectrum) = (spectruminfo(X).k0 + i - 1) * (framerate(X) / spectruminfo(X).n)
DSP.freq(i::AbstractArray{<:Real}, X::SampledSpectrum) = (spectruminfo(X).k0 .+ i .- 1) .* (framerate(X) / spectruminfo(X).n)

### transforms

function checkregion(region)
  region ∈ (1, (1,), 1:1) ||
    throw(ArgumentError("signals are transformed along time (dimension 1) only; use samples() for other regions"))
end

function gatherbins(X::PaddedSpectrum, bins)
  k0 = spectruminfo(X).k0
  p = getfield(X, :data)
  ndims(p) == 1 ? [p[k - k0 + 1] for k ∈ bins] : [p[k - k0 + 1, j] for k ∈ bins, j ∈ axes(p, 2)]
end

function incomplete(X)
  m = spectruminfo(X)
  !m.onesided && size(X, 1) == m.n &&
    throw(ArgumentError("spectrum is shifted; use ifftshift() or padded() before inverting"))
  throw(ArgumentError("spectrum is incomplete; use padded() to zero-fill it before inverting"))
end

function twosidedbins(X::SampledSpectrum)
  m = spectruminfo(X)
  m.onesided && throw(ArgumentError("spectrum is one-sided; use irfft() or brfft(), or padded() to invert"))
  size(X, 1) == m.n && mod(firstbin(X), m.n) == 0 || incomplete(X)
  samples(X)
end

twosidedbins(X::PaddedSpectrum) = gatherbins(X, 0:spectruminfo(X).n-1)

function onesidedbins(X::SampledSpectrum)
  m = spectruminfo(X)
  nbins = m.n ÷ 2 + 1
  m.onesided || return collect(selectdim(twosidedbins(X), 1, 1:nbins))
  size(X, 1) == nbins && firstbin(X) == 0 || incomplete(X)
  samples(X)
end

onesidedbins(X::PaddedSpectrum) = gatherbins(X, 0:spectruminfo(X).n÷2)

"""
    fft(x::SampledSignal)
    rfft(x::SampledSignal)

Computes the FFT of each channel of signal `x`, along time. The result is a
spectrum that keeps the sampling rate of `x`, so that [`domain`](@ref) gives the
frequency of each bin, and `ifft` (or `irfft`) restores a signal at the original
sampling rate. A spectrum may be shifted using `fftshift` and `ifftshift`, sliced,
concatenated and [`padded`](@ref).

Inverse transforms need a complete spectrum: all bins of an unshifted grid, or a
padded spectrum, from which one period is selected. `irfft` uses the non-negative
frequency bins, so it also accepts two-sided spectra.
"""
AbstractFFTs.fft(x::SampledSignal, region=1; kwargs...) =
  (checkregion(region); fftspectrum(fft(samples(x), 1; kwargs...), x))

AbstractFFTs.rfft(x::SampledSignal, region=1; kwargs...) =
  (checkregion(region); fftspectrum(rfft(samples(x), 1; kwargs...), x; onesided=true))

for f ∈ (:ifft, :bfft)
  @eval AbstractFFTs.$f(x::SampledSignal, args...; kwargs...) = $f(samples(x), args...; kwargs...)
  @eval function AbstractFFTs.$f(X::SampledSpectrum, region=1; kwargs...)
    checkregion(region)
    signal($f(twosidedbins(X), 1; kwargs...), framerate(X))
  end
end

for f ∈ (:irfft, :brfft)
  @eval AbstractFFTs.$f(x::SampledSignal, d::Integer, args...; kwargs...) = $f(samples(x), d, args...; kwargs...)
  @eval AbstractFFTs.$f(X::SampledSpectrum; kwargs...) = signal($f(onesidedbins(X), spectruminfo(X).n, 1; kwargs...), framerate(X))
  @eval function AbstractFFTs.$f(X::SampledSpectrum, d::Integer, region=1; kwargs...)
    checkregion(region)
    d == spectruminfo(X).n ||
      throw(ArgumentError("spectrum is from a $(spectruminfo(X).n)-point FFT, so cannot be inverted to $d points"))
    $f(X; kwargs...)
  end
end

for (f, from, to) ∈ ((:fftshift, :(0), :(-(m.n ÷ 2))), (:ifftshift, :(-(m.n ÷ 2)), :(0)))
  @eval function AbstractFFTs.$f(X::SampledSpectrum, dims=1)
    checkregion(dims)
    m = spectruminfo(X)
    X isa PaddedSpectrum && throw(ArgumentError("padded spectra need not be shifted, as inverse transforms select one period"))
    m.onesided && throw(ArgumentError("one-sided spectra cannot be shifted"))
    size(X, 1) == m.n || throw(ArgumentError("only complete spectra can be shifted"))
    mod(firstbin(X) - $from, m.n) == 0 ||
      throw(ArgumentError($(f == :fftshift ? "spectrum is already shifted" : "spectrum is not shifted")))
    withinfo($f(samples(X), 1), SpectrumInfo(m.fs, m.n, false, $to))
  end
end
