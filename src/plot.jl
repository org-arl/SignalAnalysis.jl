export psd, psd!, specgram

"Plot power spectral density of the signal."
function psd(s; fs=deffs, nfft=1024, plot=plot, window=nothing, label=nothing, kwargs...)
    nfft = nextfastfft(nfft)
    while nfft > size(s,1)
        nfft = div(nfft, 2)
    end
    p = welch_pgram(s, nfft; fs=fs, window=window)
    f = freq(p)
    funit = "Hz"
    if maximum(f) >= 10000
        f /= 1000.0
        funit = "kHz"
    end
    plot(f, power(p); label=label, kwargs...)
    xlabel!("Frequency ("*funit*")")
    ylabel!("Power spectral density (dB/Hz)")
end

"Plot power spectral density of the signal."
psd!(s; fs=deffs, nfft=1024, window=nothing, label=nothing, kwargs...) = psd(s; fs=fs, nfft=nfft, plot=plot!, window=window, label=label, kwargs...)

"Plot spectrogram of the signal."
function specgram(s; fs=deffs, nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, kwargs...)
    p = spectrogram(s, nfft, noverlap; fs=fs, window=window)
    t = time(p)
    f = freq(p)
    z = power(p)
    if isanalytic(s)
        ndx = sortperm(f)
        f = f[ndx]
        z = z[ndx, :]
    end
    tunit = "s"
    funit = "Hz"
    if maximum(f) >= 10000
        f /= 1000.0
        funit = "kHz"
    end
    if maximum(t) <= 1.0
        t *= 1000.0
        tunit = "ms"
    end
    heatmap(t, f, z; kwargs...)
    xlabel!("Time ("*tunit*")")
    ylabel!("Frequency ("*funit*")")
end
