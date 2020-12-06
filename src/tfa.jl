using DSP, DSP.Periodograms
using FFTW

export tfd, Wigner, Spectrogram

abstract type TFKernel end

Base.@kwdef struct Spectrogram <: TFKernel
  nfft = 256
  noverlap = 0
  window = nothing
end

Base.@kwdef struct Wigner <: TFKernel
  nfft = 0
  window = nothing
  smooth = 0
  method = :CMP2005     # alternatives: CM1980
end

struct TFD{T} <: Periodograms.TFR{T}
  power::Matrix{T}
  freq::AbstractVector
  time::AbstractVector
end

DSP.Periodograms.time(tfr::TFD) = tfr.time

"""
$(TYPEDSIGNATURES)
Computes a spectrogram time-frequency distribution of signal `s` sampled as sampling rate `fs`.

# Examples:
```julia-repl
julia> x = real.(chirp(1kHz, 10kHz, 1s, 44.1kHz));
julia> y = tfd(x, Spectrogram());
julia> plot(y)
julia> y = tfd(x, Spectrogram(nfft=512, noverlap=256, window=hamming));
julia> plot(y)
```
"""
function tfd(s, kernel::Spectrogram; onesided=eltype(s)<:Real, fs=framerate(s))
  onesided && ((eltype(s) <: Real) || (s = real.(s)))
  tfr = spectrogram(samples(s), kernel.nfft, kernel.noverlap; onesided=onesided, fs=fs, window=kernel.window)
  w = power(tfr)
  f = freq(tfr)
  if !onesided
    f = fftshift(f)
    w = fftshift(w, 1)
  end
  TFD(w, f, time(tfr))
end

"""
$(TYPEDSIGNATURES)
Computes a Wigner-Ville time-frequency distribution of signal `s` sampled as sampling rate `fs`.

# Examples:
```julia-repl
julia> x = real.(chirp(1kHz, 10kHz, 0.01s, 44.1kHz));
julia> y = tfd(x, Wigner());
julia> plot(y; clim=(0,20))
julia> y = tfd(x, Wigner(nfft=512, smooth=10, method=:CM1980, window=hamming));
julia> plot(y; clim=(0,20))
```
"""
function tfd(s, kernel::Wigner; onesided=eltype(s)<:Real, fs=framerate(s))
  N = nframes(s)
  T = eltype(s)
  nfft = kernel.nfft <= 0 ? 2N : min(kernel.nfft, 2N)
  win = kernel.window == nothing ? nothing : kernel.window isa Function ? kernel.window(nfft) : kernel.window
  x = samples(s)
  if kernel.method === :CMP2005
    # Chassande-Mottin & Pai, IEEE Sig. Proc. Letters, 12(7), 2005.
    x̂ = repeat(s; inner=2)
  elseif kernel.method === :CM1980
    # Claasen & Mecklenbräuker, Philips J. Res., 35(4/5), p276–300, 1980.
    x̂ = resample(s, 2)
    length(x̂) < 2N && (x̂ = vcat(x̂, zeros(2N-length(x̂))))
  end
  X = zeros(T, nfft, N)
  @views for n ∈ 0:N-1
    kn = min(2n, 2N-1-2n, nfft÷2-1)
    pad = nfft - length(-kn:kn)
    if pad > 0
      xx = X[floor(Int, pad/2)+2 : end-(ceil(Int, pad/2)-1), n+1] # view
    else
      xx = X[:, n+1] # view 
    end
    xx .= x̂[1 .+ 2n .+ (-kn:kn)] .* conj.(x̂[1 .+ 2n .- (-kn:kn)])
  end
  if kernel.smooth > 0
    X = copy(filt(ones(kernel.smooth)./kernel.smooth, copy(X'))') # The copying helps with memory access performance
  end
  win === nothing || (X .*= win)
  w = Array{T,2}(undef, onesided ? nfft÷2+1 : nfft, N)
  P = onesided ? plan_rfft(zeros(T, nfft)) : plan_fft(zeros(complex(T), nfft))
  wo = Array{complex(T),1}(undef, onesided ? nfft÷2+1 : nfft)
  @views for j ∈ 1:N
    mul!(wo, P, X[:,j])
    w[:,j] .= real.(wo) ./ √(2π)
  end
  if onesided
    f = FFTW.rfftfreq(nfft, fs)
  else
    f = fftshift(FFTW.fftfreq(nfft, fs))
    w = fftshift(w, 1)
  end
  TFD(w, f, (0:N-1)./fs)
end
