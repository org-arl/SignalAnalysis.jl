using FFTW: AbstractFFTs

struct SpectrumInfo
  fs::Float32
  n::Int
  onesided::Bool
  shifted::Bool
end

Base.keys(x::SpectrumInfo) = fieldnames(typeof(x))

const SampledSpectrum = MetaArray{<:Any,SpectrumInfo,T} where T

spectrum(X, m::SpectrumInfo) = MetaArray(m, X)
spectrum(X, x::SampledSignal; onesided=false) =
  spectrum(X, SpectrumInfo(framerate(x), nframes(x), onesided, false))

spectruminfo(X::SampledSpectrum) = getfield(X, :meta)

SignalBase.framerate(X::SampledSpectrum) = X.fs

samples(X::SampledSpectrum) = samples(getfield(X, :data))

function Base.show(io::IO, mime::MIME"text/plain", X::SampledSpectrum)
  m = spectruminfo(X)
  print(io, "SampledSpectrum of ", m.n, "-point signal @ ", m.fs, " Hz")
  m.onesided && print(io, ", one-sided")
  m.shifted && print(io, ", shifted")
  print(io, ", ")
  show(io, mime, samples(X))
end

"""
    domain(X::SampledSpectrum)

Returns the frequencies (in Hz) of the bins of spectrum `X`.
"""
function domain(X::SampledSpectrum)
  m = spectruminfo(X)
  k = nframes(X)
  Δf = Float64(m.fs) / m.n
  m.onesided && return (0:k-1) .* Δf
  f = fftfreq(k, k * Δf)
  m.shifted ? fftshift(f) : f
end

function checkregion(region)
  region ∈ (1, (1,), 1:1) ||
    throw(ArgumentError("signals are transformed along time (dimension 1) only; use samples() for other regions"))
end

function inverserate(X::SampledSpectrum, n)
  m = spectruminfo(X)
  m.shifted && throw(ArgumentError("spectrum is shifted; use ifftshift() before inverting"))
  n == m.n ? m.fs : Float64(m.fs) * n / m.n
end

function asserttwosided(X::SampledSpectrum)
  spectruminfo(X).onesided && throw(ArgumentError("spectrum is one-sided; use irfft() or brfft() to invert"))
  X
end

function assertonesided(X::SampledSpectrum)
  spectruminfo(X).onesided || throw(ArgumentError("spectrum is two-sided; use ifft() or bfft() to invert"))
  X
end

function onesidedlength(X::SampledSpectrum)
  n = spectruminfo(X).n
  nframes(X) == n ÷ 2 + 1 ||
    throw(ArgumentError("spectrum has been resized; specify the length of the inverse transform"))
  n
end

"""
    fft(x::SampledSignal)
    rfft(x::SampledSignal)

Computes the FFT of each channel of signal `x`, along time. The result is a
spectrum that keeps the sampling rate of `x`, so that [`domain`](@ref) gives the
frequency of each bin, and `ifft` (or `irfft`) restores a signal at the original
sampling rate. A spectrum may be shifted using `fftshift` and `ifftshift`.
"""
AbstractFFTs.fft(x::SampledSignal, region=1; kwargs...) =
  (checkregion(region); spectrum(fft(samples(x), 1; kwargs...), x))

AbstractFFTs.rfft(x::SampledSignal, region=1; kwargs...) =
  (checkregion(region); spectrum(rfft(samples(x), 1; kwargs...), x; onesided=true))

for f ∈ (:ifft, :bfft)
  @eval AbstractFFTs.$f(x::SampledSignal, args...; kwargs...) = $f(samples(x), args...; kwargs...)
  @eval function AbstractFFTs.$f(X::SampledSpectrum, region=1; kwargs...)
    checkregion(region)
    signal($f(samples(asserttwosided(X)), 1; kwargs...), inverserate(X, nframes(X)))
  end
end

for f ∈ (:irfft, :brfft)
  @eval AbstractFFTs.$f(x::SampledSignal, d::Integer, args...; kwargs...) = $f(samples(x), d, args...; kwargs...)
  @eval AbstractFFTs.$f(X::SampledSpectrum; kwargs...) = $f(X, onesidedlength(assertonesided(X)); kwargs...)
  @eval function AbstractFFTs.$f(X::SampledSpectrum, d::Integer, region=1; kwargs...)
    checkregion(region)
    signal($f(samples(assertonesided(X)), d, 1; kwargs...), inverserate(X, d))
  end
end

for (f, shifted) ∈ ((:fftshift, true), (:ifftshift, false))
  @eval function AbstractFFTs.$f(X::SampledSpectrum, dims=1)
    checkregion(dims)
    m = spectruminfo(X)
    m.onesided && throw(ArgumentError("one-sided spectra cannot be shifted"))
    m.shifted == $shifted && throw(ArgumentError($(shifted ? "spectrum is already shifted" : "spectrum is not shifted")))
    spectrum($f(samples(X), 1), SpectrumInfo(m.fs, m.n, m.onesided, $shifted))
  end
end

"""
    fftfreq(x::SampledSignal)
    fftfreq(n, x::SampledSignal)

Frequencies of the bins of an `n`-point FFT of signal `x`. If `n` is not
specified, the number of frames in `x` is used.
"""
AbstractFFTs.fftfreq(x::SampledSignal) = fftfreq(nframes(x), framerate(x))
AbstractFFTs.fftfreq(n::Integer, x::SampledSignal) = fftfreq(n, framerate(x))
