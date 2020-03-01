export samplingrate, isanalytic, analytic
export energy, meantime, rmsduration, meanfrequency, bandwidth, ifreq

deffs = 1.0

"Set default sampling rate."
function samplingrate(fs)
  global deffs
  deffs = freqQ(fs)
end

"Get default sampling rate."
samplingrate() = deffs

"Convert a signal to analytic representation."
analytic(s::AbstractArray) = isanalytic(s) ? s : hilbert(s)

"Check if signal is analytic."
isanalytic(s::AbstractArray) = eltype(s) <: Complex

"Get time vector corresponding to each sample in signal."
Base.time(s::AbstractArray; fs=deffs) = float(0:size(s,1)-1)./freqQ(fs)

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

"Get bandwidth of the signal."
function bandwidth(s; fs=deffs, nfft=1024, window=nothing)
  mapslices(s; dims=1) do s1
    p = welch_pgram(s1, ceil(Int, length(s1)/nfft); fs=freqQ(fs), window=window)
    f = freq(p)
    f0 = wmean(f, power(p))
    sqrt.(wmean((f.-f0).^2, power(p)))
  end
end

### utility functions

wmean(x, w::AbstractVector) = sum(x.*w) ./ sum(w)
wmean(x, w::AbstractMatrix) = sum(x.*w; dims=1) ./ sum(w; dims=1)
