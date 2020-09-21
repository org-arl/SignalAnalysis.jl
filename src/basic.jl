export energy, meantime, rmsduration, meanfrequency, rmsbandwidth, ifrequency

"""
$(SIGNATURES)
Computes total signal energy.
"""
energy(s::AbstractVector; fs=framerate(s)) = sum(abs2, s)/inHz(fs)
energy(s::AbstractMatrix; fs=framerate(s)) = vec(sum(abs2, s; dims=1))./inHz(fs)

"""
$(SIGNATURES)
Computes mean time of the signal.
"""
meantime(s::SampledSignal) = wmean(domain(s), abs2.(samples(s)))
meantime(s; fs) = wmean((0:size(s,1)-1)/fs, abs2.(s))

"""
$(SIGNATURES)
Computes RMS duration of the signal.
"""
rmsduration(s::SampledSignal) = sqrt.(wmean(domain(s).^2, abs2.(samples(s))) .- meantime(s).^2)
rmsduration(s; fs) = sqrt.(wmean(((0:size(s,1)-1)/fs).^2, abs2.(s)) .- meantime(s; fs=fs).^2)

"""
$(SIGNATURES)
Computes instantaneous frequency of the signal.
"""
function ifrequency(s; fs=framerate(s))
  s1 = analytic(s)
  f1 = inHz(fs)/(2π) * diff(unwrap(angle.(s1); dims=1); dims=1)
  f2 = Array{Float32}(undef, size(s))
  f2[1,:] = f1[1,:]
  f2[2:end-1,:] .= (f1[1:end-1,:] .+ f1[2:end,:])./2
  f2[end,:] = f1[end,:]
  signal(f2, fs)
end

"""
$(SIGNATURES)
Computes mean frequency of a signal.
"""
function meanfrequency(s::AbstractMatrix; fs=framerate(s), nfft=1024, window=nothing)
  fs = inHz(fs)
  rv = mapslices(samples(s); dims=1) do s1
    p = welch_pgram(s1, min(nfft, length(s1)); fs=fs, window=window)
    f = freq(p)
    wmean(f, power(p))
  end
  dropdims(rv, dims=1)
end

function meanfrequency(s::AbstractVector; fs=framerate(s), nfft=1024, window=nothing)
  p = welch_pgram(s, min(nfft, length(s)); fs=inHz(fs), window=window)
  f = freq(p)
  wmean(f, power(p))
end

"""
$(SIGNATURES)
Computes RMS bandwidth of a signal.
"""
function rmsbandwidth(s::AbstractMatrix; fs=framerate(s), nfft=1024, window=nothing) where T
  fs = inHz(fs)
  rv = mapslices(samples(s); dims=1) do s1
    p = welch_pgram(s1, min(nfft, length(s1)); fs=fs, window=window)
    f = freq(p)
    f0 = wmean(f, power(p))
    sqrt.(wmean((f.-f0).^2, power(p)))
  end
  dropdims(rv, dims=1)
end

function rmsbandwidth(s::AbstractVector; fs=framerate(s), nfft=1024, window=nothing) where T
  p = welch_pgram(s, min(nfft, length(s)); fs=inHz(fs), window=window)
  f = freq(p)
  f0 = wmean(f, power(p))
  sqrt.(wmean((f.-f0).^2, power(p)))
end

### utility functions

wmean(x::AbstractVector, w::AbstractVector) = (x'w) / sum(w)
wmean(x, w::AbstractVector) = sum(x.*w) ./ sum(w)
wmean(x, w::AbstractMatrix) = dropdims(sum(x.*w; dims=1) ./ sum(w; dims=1), dims=1)
