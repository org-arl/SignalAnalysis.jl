export tfd, Wigner, Spectrogram

using DSP.Periodograms
using FFTW

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

function tfd(s, kernel::Wigner; onesided=eltype(s)<:Real, fs=framerate(s))
  N = nframes(s)
  nfft = kernel.nfft <= 0 ? 2N : min(kernel.nfft, 2N)
  win = kernel.window == nothing ? nothing : kernel.window isa Function ? kernel.window(nfft) : kernel.window
  x = samples(s)
  if kernel.method == :CMP2005
    # Chassande-Mottin & Pai, IEEE Sig. Proc. Letters, 12(7), 2005.
    x̂ = repeat(s; inner=2)
  elseif kernel.method == :CM1980
    # Claasen & Mecklenbräuker, Philips J. Res., 35(4/5), p276–300, 1980.
    x̂ = resample(s, 2)
    length(x̂) < 2N && (x̂ = vcat(x̂, zeros(2N-length(x̂))))
  end
  X = Array{Float64, 2}(undef, N, nfft)
  for n ∈ 0:N-1
    kn = min(2n, 2N-1-2n, nfft÷2-1)
    xx = x̂[1 .+ 2n .+ (-kn:kn)] .* conj(x̂[1 .+ 2n .- (-kn:kn)])
    pad = nfft - length(xx)
    pad > 0 && (xx = vcat(zeros(floor(Int, pad/2)+1), xx, zeros(ceil(Int, pad/2)-1)))
    onesided && ((eltype(s) <: Real) || (xx = real.(xx)))
    X[n+1,:] .= xx
  end
  kernel.smooth > 0 && (X = filt(ones(kernel.smooth)./kernel.smooth, X))
  win == nothing || (X = X .* win')
  w = Array{Float64,2}(undef, onesided ? nfft÷2+1 : nfft, N)
  P = onesided ? plan_rfft(zeros(Float64, nfft)) : plan_fft(zeros(ComplexF64, nfft))
  for j ∈ 1:N
    w[:,j] = real.(P * X[j,:]) / √(2π)
  end
  if onesided
    f = FFTW.rfftfreq(nfft, fs)
  else
    f = fftshift(FFTW.fftfreq(nfft, fs))
    w = fftshift(w, 1)
  end
  TFD(w, f, (0:N-1)./fs)
end
