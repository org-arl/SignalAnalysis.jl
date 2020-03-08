export samplingrate, isanalytic, analytic
export padded, slide
export energy, meantime, rmsduration, meanfrequency, rmsbandwidth, ifreq, fir

deffs = 1.0

"Set default sampling rate. Returns previous sampling rate."
function samplingrate(fs)
  global deffs
  old = deffs
  deffs = freqQ(fs)
  return old
end

"Get default sampling rate."
samplingrate() = deffs

"Run code block with temporary sampling rate."
function samplingrate(f::Function, fs)
  old = samplingrate(fs)
  rv = f()
  samplingrate(old)
  return rv
end

"Generate a padded view of a signal with optional delay/advance."
function padded(s::AbstractVector{T}, padding; delay=0, fill=zero(T)) where {T, N}
  if length(padding) == 1
    left = padding
    right = padding
  else
    left = padding[1]
    right = padding[2]
  end
  PaddedView(fill, s, (1-left:length(s)+right,), (1+delay:delay+length(s),))
end

"Slide a window over a signal, process each window."
function slide(f::Function, s::AbstractVector, nsamples, overlap=0, args...; showprogress=true)
  @assert overlap < nsamples "overlap must be less than nsamples"
  n = size(s,1)
  m = nsamples - overlap
  mmax = (n-nsamples)÷m
  showprogress && (p = Progress(mmax+1, 1, "Processing: "))
  for j = 0:mmax
    s1 = @view s[j*m+1:j*m+nsamples]
    f(s1, j+1, j*m+1, args...)
    showprogress && next!(p)
  end
end

"Slide a window over a signal, process each window, and collect the results."
function slide(f::Function, ::Type{T}, s::AbstractVector, nsamples, overlap=0, args...; showprogress=true) where {T}
  @assert overlap < nsamples "overlap must be less than nsamples"
  n = size(s,1)
  m = nsamples - overlap
  mmax = (n-nsamples)÷m
  out = Array{T,1}(undef, 1+mmax)
  showprogress && (p = Progress(mmax+1, 1, "Processing: "))
  for j = 0:mmax
    s1 = @view s[j*m+1:j*m+nsamples]
    out[j+1] = f(s1, j+1, j*m+1, args...)
    showprogress && next!(p)
  end
  return out
end

"Convert a signal to analytic representation."
analytic(s::AbstractArray) = isanalytic(s) ? s : hilbert(s)

"Check if signal is analytic."
isanalytic(s::AbstractArray) = eltype(s) <: Complex

"Get time vector corresponding to each sample in signal."
Base.time(s::AbstractArray; fs=deffs, t0=0.0) = t0 .+ float(0:size(s,1)-1)./freqQ(fs)

"Convert time to index."
toindex(t; fs=deffs) = 1 + round(Int, timeQ(t)*freqQ(fs))

"Get total signal energy."
energy(s::AbstractVector; fs=deffs) = sum(abs2, s)/freqQ(fs)
energy(s::AbstractMatrix; fs=deffs) = dropdims(sum(abs2, s; dims=1); dims=1)./freqQ(fs)

"Get mean time of the signal."
meantime(s; fs=deffs) = wmean(time(s; fs=fs), abs.(s).^2)

"Get RMS duration of the signal."
rmsduration(s; fs=deffs) = sqrt.(wmean(time(s; fs=fs).^2, abs.(s).^2) .- meantime(s; fs=fs).^2)

"Get instantaneous frequency of the signal."
function ifreq(s; fs=deffs)
  s1 = analytic(s)
  f1 = freqQ(fs)/(2π) * diff(unwrap(angle.(s1); dims=1); dims=1)
  vcat(f1[[1],:], (f1[1:end-1,:]+f1[2:end,:])/2, f1[[end],:])
end

"Get mean frequency of the signal."
function meanfrequency(s; fs=deffs, nfft=1024, window=nothing)
  mapslices(s; dims=1) do s1
    p = welch_pgram(s1, ceil(Int, length(s1)/nfft); fs=freqQ(fs), window=window)
    f = freq(p)
    wmean(f, power(p))
  end
end

"Get RMS bandwidth of the signal."
function rmsbandwidth(s; fs=deffs, nfft=1024, window=nothing)
  mapslices(s; dims=1) do s1
    p = welch_pgram(s1, ceil(Int, length(s1)/nfft); fs=freqQ(fs), window=window)
    f = freq(p)
    f0 = wmean(f, power(p))
    sqrt.(wmean((f.-f0).^2, power(p)))
  end
end

"Design FIR filter."
function fir(n, f1, f2=nothing; fs=deffs, method=FIRWindow(hanning(n)))
  fs = freqQ(fs)
  if f1 == 0
    f = Lowpass(freqQ(f2); fs=fs)
  elseif f2 == nothing || freqQ(f2) == fs/2
    f = Highpass(freqQ(f1); fs=fs)
  else
    f = Bandpass(freqQ(f1), freqQ(f2); fs=fs)
  end
  return digitalfilter(f, method)
end

### utility functions

wmean(x, w::AbstractVector) = sum(x.*w) ./ sum(w)
wmean(x, w::AbstractMatrix) = sum(x.*w; dims=1) ./ sum(w; dims=1)
