using FFTW: AbstractFFTs
import MetaArrays

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

Base.getindex(p::PeriodicView{T,1}, i::Int) where T =
  checkbounds(Bool, p.data, i) ? @inbounds(p.data[i]) : periodicvalue(p, i)
Base.getindex(p::PeriodicView{T,2}, i::Int, j::Int) where T =
  checkbounds(Bool, p.data, i, j) ? @inbounds(p.data[i, j]) : periodicvalue(p, i, j)

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

function sidedness(X, m::SpectrumInfo, i)
  X isa PaddedSpectrum && m.onesided || return m
  len = i isa Colon ? size(X, 1) : length(i)
  0 ≤ m.k0 && m.k0 + len - 1 ≤ m.n ÷ 2 && return m
  SpectrumInfo(m.fs, m.n, false, m.k0)
end

function spectrumindex(f, X, I...)
  J = map(i -> toindex(X, i), I)
  data = getfield(X, :data)
  m = sliceinfo(X, J...)
  m === nothing ? f(data, J...) : withinfo(f(data, J...), sidedness(X, m, first(J)))
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
domain(X::SampledSpectrum) = (firstbin(X):lastbin(X)) .* (Float64(framerate(X)) / spectruminfo(X).n)

"""
    fftfreq(X::SampledSpectrum)

Returns the frequencies (in Hz) of the bins of spectrum `X`, wrapped to the
interval `[-fs/2, fs/2)` for two-sided spectra.
"""
function AbstractFFTs.fftfreq(X::SampledSpectrum)
  f = collect(domain(X))
  spectruminfo(X).onesided && return f
  fs = Float64(framerate(X))
  mod.(f .+ fs / 2, fs) .- fs / 2
end

"""
    fftfreq(x::SampledSignal)
    fftfreq(n, x::SampledSignal)

Frequencies of the bins of an `n`-point FFT of signal `x`. If `n` is not
specified, the number of frames in `x` is used.
"""
AbstractFFTs.fftfreq(x::SampledSignal) = fftfreq(nframes(x), Float64(framerate(x)))
AbstractFFTs.fftfreq(n::Integer, x::SampledSignal) = fftfreq(n, Float64(framerate(x)))

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
DSP.freq(i::Real, X::SampledSpectrum) = (spectruminfo(X).k0 + i - 1) * (Float64(framerate(X)) / spectruminfo(X).n)
DSP.freq(i::AbstractArray{<:Real}, X::SampledSpectrum) = (spectruminfo(X).k0 .+ i .- 1) .* (Float64(framerate(X)) / spectruminfo(X).n)

### values between bins

const FrequencyLike = Union{Real,Units.Unitful.Frequency}

"""
    (X::SampledSpectrum)(f)
    (X::SampledSpectrum)(f, ch)

Evaluate spectrum `X` at frequency `f`, which may lie between bins. The frequency
`f` is in Hz, or may be given with frequency units (e.g. `1.5kHz`). The result is
the discrete-time Fourier transform of the signal that the spectrum represents, so
it is exact, and agrees with the bins at bin frequencies, up to rounding. As for
the bins, time 0 of that signal is its first sample. The result is periodic in
frequency, with period `fs`, and one-sided spectra are extended by conjugate
symmetry, as for a real signal.

The spectrum must be complete: all bins of a two-sided spectrum (shifted or not),
all non-negative frequency bins of a one-sided spectrum, or a [`padded`](@ref)
spectrum, whose missing bins are zero. Each value takes `O(n)` operations for an
`n`-point FFT, so for values on a fine regular grid, zero-padding the signal
before the FFT is faster.

For a single-channel spectrum, `X(f)` returns a single value. For a multichannel
spectrum, it returns a vector with one value per channel, and `X(f, ch)` returns
the value of channel `ch`. If `f` is a vector of frequencies, the result is a
vector (for a single-channel spectrum, or when `ch` is given), or a matrix with
one row per frequency and one column per channel.

# Examples:
```julia-repl
julia> X = fft(signal(cos.(2π * 1050 * (0:99) / 8000), 8kHz));

julia> abs(X(1kHz))
22.887595468624955

julia> abs(X(1050))
50.48110116186148
```
"""
function (X::SampledSpectrum)(f::FrequencyLike)
  B, ks = binswithindices(X)
  dtftvalue(B, ks, fractionalbin(f, X))
end

function (X::SampledSpectrum)(f::FrequencyLike, ch::Integer)
  Base.checkindex(Bool, axes(X, 2), ch) || throw(BoundsError(X, (f, ch)))
  B, ks = binswithindices(X)
  dtftvalue(B, ks, fractionalbin(f, X), ch)
end

function (X::SampledSpectrum)(f::AbstractVector)
  B, ks = binswithindices(X)
  dtftvalue(B, ks, fractionalbin.(f, Ref(X)))
end

function (X::SampledSpectrum)(f::AbstractVector, ch::Integer)
  Base.checkindex(Bool, axes(X, 2), ch) || throw(BoundsError(X, (f, ch)))
  B, ks = binswithindices(X)
  [dtftvalue(B, ks, fractionalbin(f1, X), ch) for f1 ∈ f]
end

fractionalbin(f, X) = Float64(inHz(f)) * spectruminfo(X).n / Float64(framerate(X))

# bins of a complete spectrum, with the FFT bin number of each row
function binswithindices(X::SampledSpectrum)
  m = spectruminfo(X)
  X isa PaddedSpectrum && return gatherbins(X, 0:m.n-1), 0:m.n-1
  m.onesided && return gatherbins(padded(withinfo(onesidedbins(X), m)), 0:m.n-1), 0:m.n-1
  size(X, 1) == m.n || incomplete(X)
  samples(X), firstbin(X):lastbin(X)
end

dtftvalue(B::AbstractVector, ks, ν::Real) = dtftvalue(B, ks, ν, 1)
dtftvalue(B::AbstractMatrix, ks, ν::Real) = [dtftvalue(B, ks, ν, j) for j ∈ axes(B, 2)]
dtftvalue(B::AbstractVector, ks, νs::AbstractVector) = [dtftvalue(B, ks, ν, 1) for ν ∈ νs]
dtftvalue(B::AbstractMatrix, ks, νs::AbstractVector) = [dtftvalue(B, ks, ν, j) for ν ∈ νs, j ∈ axes(B, 2)]

function dtftvalue(B, ks, ν::Real, ch)
  T = complex(float(eltype(B)))
  n = length(ks)
  acc = zero(T)
  for (i, k) ∈ zip(axes(B, 1), ks)
    acc += T(dirichlet(ν - k, n)) * B[i, ch]
  end
  acc
end

# DTFT at fractional bin δ of a unit bin at 0 of an n-point DFT
function dirichlet(δ, n)
  isinteger(δ / n) && return one(ComplexF64)
  cispi(-δ * (n - 1) / n) * sinpi(δ) / (n * sinpi(δ / n))
end

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
  throw(ArgumentError("spectrum is incomplete; use padded() to zero-fill the missing bins"))
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

### reductions

const SampledArray = Union{SampledSignal,SampledSpectrum}

Base.mapreduce(f, op, A::SampledArray; dims=:, kwargs...) = mapreduce(f, op, samples(A); dims=dims, kwargs...)
Base.findmax(f, A::SampledArray; dims=:) = findmax(f, samples(A); dims=dims)
Base.findmin(f, A::SampledArray; dims=:) = findmin(f, samples(A); dims=dims)
Base.findmax(A::SampledArray; dims=:) = findmax(samples(A); dims=dims)
Base.findmin(A::SampledArray; dims=:) = findmin(samples(A); dims=dims)
Statistics.mean(A::SampledArray; dims=:) = mean(samples(A); dims=dims)
Statistics.mean(f, A::SampledArray; dims=:) = mean(f, samples(A); dims=dims)
Statistics.var(A::SampledArray; dims=:, kwargs...) = var(samples(A); dims=dims, kwargs...)
Statistics.std(A::SampledArray; dims=:, kwargs...) = std(samples(A); dims=dims, kwargs...)
Statistics.median(A::SampledArray; dims=:) = median(samples(A); dims=dims)

### combining

function MetaArrays.metamerge(a::SpectrumInfo, b::SpectrumInfo)
  a == b && return a
  a.fs == b.fs || throw(ArgumentError("cannot combine spectra with different sampling rates ($(a.fs) and $(b.fs) Hz)"))
  a.n == b.n || throw(ArgumentError("cannot combine spectra from FFTs of different lengths ($(a.n) and $(b.n) points)"))
  a.onesided == b.onesided || throw(ArgumentError("cannot combine a one-sided spectrum with a two-sided one"))
  throw(ArgumentError("cannot combine spectra covering different bins (starting at bins $(a.k0) and $(b.k0)); use samples() on one of them to keep the bins of the other"))
end

function MetaArrays.metamerge(a::SamplingInfo, b::SamplingInfo)
  a == b && return a
  throw(ArgumentError("cannot combine signals with different sampling rates ($(a.fs) and $(b.fs) Hz)"))
end

MetaArrays.metamerge(::SpectrumInfo, ::SamplingInfo) =
  throw(ArgumentError("cannot combine a spectrum with a signal; use samples() if this is intended"))
MetaArrays.metamerge(::SamplingInfo, ::SpectrumInfo) =
  throw(ArgumentError("cannot combine a signal with a spectrum; use samples() if this is intended"))

Base.diff(X::SampledSpectrum; dims::Integer) = diff(samples(X); dims=dims)
Base.diff(X::MetaArray{<:AbstractVector,SpectrumInfo}) = diff(samples(X))
