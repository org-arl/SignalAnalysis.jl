export energy, meantime, rmsduration, meanfrequency, bandwidth, ifreq, psd, specgram

"Get total signal energy."
energy(s) = energy(sampledsignal(s))
energy(s::AxisArray{<:Any,1,<:Any}) = sum(abs2, s) / samplingrate(s)
energy(s::AxisArray{<:Any,2,<:Any}) = sum(abs2, s; dims=1).data ./ samplingrate(s)

"Get mean time of the signal."
meantime(s) = meantime(sampledsignal(s))
meantime(s::AxisArray) = wmean(time(s), abs.(s).^2)

"Get RMS duration of the signal."
rmsduration(s) = rmsduration(sampledsignal(s))
rmsduration(s::AxisArray) = sqrt.(wmean(time(s).^2, abs.(s).^2) .- meantime(s).^2)

"Get instantaneous frequency of the signal."
function ifreq(s)
    s1 = sampledsignal(s, analytic=true)
    f1 = samplingrate(s1)/(2π) * diff(unwrap(angle.(s1); dims=1); dims=1)
    f = vcat(f1[[1],:], (f1[1:end-1,:]+f1[2:end,:])/2, f1[[end],:])
    AxisArray(f, s1.axes...)
end

"Get power spectral density of the signal."
psd(s; nfft=1024, kwargs...) = psd(sampledsignal(s), nfft=nfft, kwargs...)
function psd(s::AxisArray; nfft=1024, kwargs...)
    fs = samplingrate(s)
    nfft = nextfastfft(nfft)
    while nfft > size(s,1)
        nfft = div(nfft, 2)
    end
    p = welch_pgram(s, nfft, fs=fs, kwargs...)
    AxisArray(power(p), Axis{:frequency}(freq(p)))
end

"Get spectrogram of the signal."
specgram(s; nfft=min(div(length(s),8),256), noverlap=div(n,2), kwargs...) = specgram(sampledsignal(s), nfft, noverlap, kwargs...)
function specgram(s::AxisArray; nfft=min(div(length(s),8),256), noverlap=div(n,2), kwargs...)
    fs = samplingrate(s)
    # FIXME: spectrogram does not work correctly with analytic signal input, so forcing signal to be real
    p = spectrogram(sampledsignal(s, analytic=false), nfft, noverlap, fs=fs, kwargs...)
    AxisArray(power(p), Axis{:frequency}(freq(p)), Axis{:time}(time(p)))
end

"Get mean frequency of the signal."
meanfrequency(s; nfft=1024, window=nothing) = meanfrequency(signal(s); nfft=nfft, window=window)
function meanfrequency(s::AxisArray; nfft=1024, window=nothing)
    fs = samplingrate(s)
    mapslices(s; dims=1) do s1
        p = welch_pgram(s1, ceil(Int, length(s1)/nfft), fs=fs, window=window)
        f = freq(p)
        wmean(f, power(p))
    end
end

"Get bandwidth of the signal."
bandwidth(s; nfft=1024, window=nothing) = bandwidth(signal(s); nfft=nfft, window=window)
function bandwidth(s::AxisArray; nfft=1024, window=nothing)
    fs = samplingrate(s)
    mapslices(s; dims=1) do s1
        p = welch_pgram(s1, ceil(Int, length(s1)/nfft), fs=fs, window=window)
        f = freq(p)
        f0 = wmean(f, power(p))
        sqrt.(wmean((f.-f0).^2, power(p)))
    end
end

"Get signal envelope."
Base.abs(s::AxisArray) = AxisArray(abs.(s.data), s.axes...)

"Get signal power."
Base.abs2(s::AxisArray) = AxisArray(abs2.(s.data), s.axes...)

### utility functions

wmean(x, w::AbstractVector) = sum(x.*w) ./ sum(w)
wmean(x, w::AbstractMatrix) = sum(x.*w; dims=1) ./ sum(w; dims=1)
